"""Interpret peptide and molecular inputs beneath display, CLI and HTTP adapters.

The result records the input format separately from normalized CABILN. In
particular, Converter.get_biln() emits CABILN with modern attachment slots; it
must never be fed back through the converter's legacy-only BILN reader.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class MolecularInput:
    format: str
    source: str = ""
    sequence: object = None
    assembly: object = None
    molecule: object = None

    def assemble(self, depiction="local"):
        """Assemble once within this operation, after any CLI PDB-name changes."""
        if self.molecule is None and self.sequence is not None:
            from pyPept.molecule import Molecule

            self.assembly = Molecule(self.sequence, depiction=depiction)
            self.molecule = self.assembly.get_molecule(fmt="ROMol")
        return self.molecule


@dataclass(frozen=True)
class FormattedSource:
    """Verified notation and its target-to-source monomer correspondence."""

    source: str
    text: str
    occurrence_order: tuple[int, ...]
    sequence: object
    peptide: object
    warnings: tuple[str, ...]


def read_input(source, *, input_format="cabiln", warning_sink=None, track_source=False):
    """Interpret an explicit format, retaining objects supplied by Python callers.

    Legacy BILN and CHUCKLES retain R3 sidechain meaning. HELM's converter output
    is already CABILN; only the original legacy input needs slot normalization.
    Parsing alone lets the CLI assign PDB names before its chosen depiction.
    """
    from rdkit import Chem

    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence
    from pyPept.notation_conversion import biln_to_cabiln

    if isinstance(source, Chem.rdchem.Mol):
        return MolecularInput("MOL", molecule=source)
    if isinstance(source, Molecule):
        return MolecularInput(
            "MOL", assembly=source, molecule=source.get_molecule(fmt="ROMol")
        )
    if isinstance(source, Sequence):
        return MolecularInput("CABILN", source.s_inputbiln, sequence=source)
    if not isinstance(source, str):
        raise TypeError(
            "Input must be notation, a Sequence, a Molecule, or an RDKit molecule"
        )

    kind = input_format.lower()
    if kind == "smiles":
        molecule = Chem.MolFromSmiles(source)
        if molecule is None:
            raise ValueError("Invalid SMILES")
        return MolecularInput("SMILES", source, molecule=molecule)
    if kind == "mol":
        molecule = Chem.MolFromMolBlock(source)
        if molecule is None:
            raise ValueError("Invalid .mol data")
        return MolecularInput("MOL", source, molecule=molecule)
    if kind == "cabiln":
        sequence, normalized = parse_source(
            source, warning_sink, track_source=track_source
        )
    else:
        if kind == "biln":
            normalized = biln_to_cabiln(source)
        elif kind in ("helm", "chuckles"):
            from pyPept.converter import Converter

            normalized = Converter(**{kind: source}).get_biln()
        elif kind == "fasta":
            normalized = "-".join(source)
        else:
            raise ValueError(f"Unknown input format: {input_format!r}")
        sequence = Sequence(
            normalized, warning_sink=warning_sink, track_source=track_source
        )
    return MolecularInput(kind.upper(), normalized, sequence=sequence)


def detect_input(
    source, *, policy="reference", on_notation_error=None, depiction="local"
):
    """Apply the established reference, conversion, or display precedence.

    Reference and conversion prefer SMILES, then HELM, then legacy BILN.
    Reference additionally tries positional CABILN. Display prefers peptide
    notation for ambiguous bare tokens, while SMILES punctuation gets first try.
    An invalid punctuated SMILES can still be valid CABILN. Callers that draw
    their own coordinates can skip the initial layout with ``depiction=None``.
    """
    if not isinstance(source, str):
        parsed = read_input(source)
        parsed.assemble(depiction=depiction)
        return parsed
    text = source if policy == "display" else source.strip()
    if policy == "display":
        formats = (
            ("smiles", "cabiln")
            if set("()[]=#@+\\/") & set(text)
            else ("cabiln", "smiles")
        )
    elif policy in ("reference", "conversion"):
        formats = ("smiles", "helm", "biln")
        if policy == "reference":
            formats += ("cabiln",)
    else:
        raise ValueError(f"Unknown input policy: {policy!r}")
    for kind in formats:
        if kind == "helm" and not ("PEPTIDE" in text.upper() and "$" in text):
            continue
        try:
            import warnings

            parsed = read_input(
                text,
                input_format=kind,
                track_source=policy == "conversion",
                warning_sink=warnings.warn if policy == "display" else None,
            )
            if policy != "conversion":
                parsed.assemble(depiction=depiction)
                if parsed.molecule is None:
                    raise ValueError("Assembly produced no molecule")
            return parsed
        except Exception as exc:
            if policy == "conversion" and kind == "helm":
                raise ValueError(f"HELM parse failed: {exc}") from exc
            if (
                policy == "display"
                and kind == formats[0] == "cabiln"
                and on_notation_error
            ):
                on_notation_error(exc)
    if policy == "conversion":
        raise ValueError("Could not parse input as SMILES, BILN, or HELM")
    raise ValueError("Could not parse as SMILES, BILN, HELM, or CABILN")


def convert_input(source, notation="percent", *, input_format="auto"):
    """Convert the requested format, or use the established auto precedence."""
    parsed = (
        detect_input(source, policy="conversion")
        if input_format.lower() == "auto"
        else read_input(source, input_format=input_format, track_source=True)
    )
    if parsed.format == "SMILES":
        from pyPept.smiles import convert_smiles

        return parsed.format, convert_smiles(parsed.source, notation=notation)
    try:
        return parsed.format, _format_parsed(parsed, notation).text
    except Exception as exc:
        if parsed.format == "HELM":
            raise ValueError(f"HELM parse failed: {exc}") from exc
        raise ValueError("Could not parse input as SMILES, BILN, or HELM") from exc


def split_top_level(text: str, delimiter: str) -> list[str]:
    """Split notation while preserving brackets and inline SMILES tokens."""
    parts, stack, start = [], [], 0
    closing = {"[": "]", "{": "}", "(": ")", "<": ">"}
    for index, char in enumerate(text):
        if stack and stack[-1] == ">":
            if char == ">":
                stack.pop()
        elif char in closing:
            stack.append(closing[char])
        elif char in closing.values():
            if not stack or stack.pop() != char:
                raise ValueError("Unbalanced brackets in CABILN")
        elif char == delimiter and not stack:
            parts.append(text[start:index])
            start = index + 1
    if stack:
        raise ValueError("Unbalanced brackets in CABILN")
    parts.append(text[start:])
    return parts


def _to_bracket(cabiln: str) -> str:
    """Normalise positional-% branch notation to bracket form for Sequence.

    Pure crosslink branches (%TBMB.!1.!2.!3, %C20FA-AEEA-E_g.!1, etc.) are
    already parseable by Sequence and are left untouched.  Only branches where
    every segment lacks an !n marker need conversion.
    """
    if "%" not in cabiln:
        return cabiln
    branch_parts = split_top_level(cabiln, "%")[1:]
    if all("!" in p for p in branch_parts):
        return cabiln  # all crosslink — Sequence handles it
    from pyPept.notation_conversion import cabiln_to_bracket

    return cabiln_to_bracket(cabiln)


def parse_source(cabiln: str, warning_sink=None, *, track_source=False):
    """Return the parsed source explicitly, including legacy positional input.

    Callers exposing residue selections must use the returned source for edits;
    legacy conversion can change instance ordering.
    """
    from pyPept.sequence import Sequence
    from pyPept.source import SourceError

    messages = []
    parsed_source = cabiln
    try:
        sequence = Sequence(
            cabiln, warning_sink=messages.append, track_source=track_source
        )
    except ValueError:
        parsed_source = _to_bracket(cabiln)
        if parsed_source == cabiln:
            raise
        messages.clear()
        try:
            sequence = Sequence(
                parsed_source, warning_sink=messages.append, track_source=track_source
            )
        except SourceError as error:
            # The fallback rewrites source positions before parsing. Its spans
            # cannot select the original input that the user still sees.
            raise SourceError(str(error), hint=error.hint) from error
    if warning_sink is not None:
        for message in messages:
            warning_sink(message)
    return sequence, parsed_source


def format_source(source: str, notation: str, *, canonical: bool = False) -> str:
    """Format resolved occurrences and verify their identity and chemistry.

    Canonical output ignores source layout and spelling, retaining the selected
    monomer decomposition. The default preserves existing layout preferences.
    """
    return format_source_details(source, notation, canonical=canonical).text


def format_source_details(
    source: str, notation: str, *, canonical: bool = False
) -> FormattedSource:
    """Retain the validated target parse for presentation without reparsing it."""
    parsed = read_input(source, track_source=True)
    return _format_parsed(parsed, notation, canonical=canonical)


def _format_parsed(parsed, notation, *, canonical=False):
    from rdkit import Chem

    from pyPept.molecule import Molecule
    from pyPept.peptide import Connection, Endpoint, Peptide, serialize
    from pyPept.sequence import Sequence
    from pyPept.structure import compare_structures

    sequence = parsed.sequence
    peptide = Peptide.from_sequence(sequence)
    emission = serialize(peptide, notation=notation, canonical=canonical)
    if (
        len(emission.occurrence_order) != len(peptide.occurrences)
        or set(emission.occurrence_order) != set(range(len(peptide.occurrences)))
    ):
        raise ValueError("Notation conversion would lose or duplicate a monomer")
    messages = []
    converted = Sequence(
        emission.text, warning_sink=messages.append, track_source=True
    )
    if len(converted.s_monomers) != len(emission.occurrence_order):
        raise ValueError("Notation conversion would lose or duplicate a monomer")
    target = Peptide.from_sequence(converted)
    remapped_connections = {
        Connection(*(
            Endpoint(emission.occurrence_order[end.occurrence_id], end.slot)
            for end in edge.endpoints
        ))
        for edge in target.connections
    }
    if set(peptide.connections) != remapped_connections:
        raise ValueError("Notation conversion would change monomer connections")
    for index, identity in enumerate(emission.occurrence_order):
        before = sequence.s_monomers[identity]
        after = converted.s_monomers[index]
        if (
            Chem.MolToSmiles(before["m_romol"]) != Chem.MolToSmiles(after["m_romol"])
            or before["m_Rgroups"] != after["m_Rgroups"]
        ):
            raise ValueError(
                "Notation conversion would change the molecular structure of a monomer"
            )
    if parsed.molecule is None:
        parsed.assembly = Molecule(peptide, depiction=None)
        parsed.molecule = parsed.assembly.get_molecule(fmt="ROMol")
    product = Molecule(target, depiction=None).get_molecule(fmt="ROMol")
    if not compare_structures(parsed.molecule, product).exact:
        raise ValueError("Notation conversion would change the molecular structure")
    return FormattedSource(
        parsed.source, emission.text, emission.occurrence_order,
        converted, target, tuple(messages),
    )
