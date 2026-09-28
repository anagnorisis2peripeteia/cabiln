"""Notation normalization shared by builder and conversion requests."""

from __future__ import annotations

import re


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


def backbone_token_indices(tokens: list[str]) -> list[int]:
    """Terminal !n markers are bonds, and do not consume a residue index."""
    return [
        i
        for i, token in enumerate(tokens)
        if token and not re.fullmatch(r"![A-Za-z0-9_]+", token)
    ]


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
    from pyPept.sequence import cabiln_to_bracket

    return cabiln_to_bracket(cabiln)


def parse_source(cabiln: str, warning_sink=None, *, track_source=False):
    """Return the parsed source explicitly, including legacy positional input.

    Callers exposing residue selections must use the returned source for edits;
    legacy conversion can change instance ordering.
    """
    from pyPept.sequence import Sequence

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
        sequence = Sequence(
            parsed_source, warning_sink=messages.append, track_source=track_source
        )
        messages.insert(0, "Converted legacy positional notation to bracket form.")
    if warning_sink is not None:
        for message in messages:
            warning_sink(message)
    return sequence, parsed_source


def format_source(source: str, notation: str) -> str:
    """Format resolved occurrences and verify their identity and chemistry."""
    from rdkit import Chem

    from pyPept.molecule import Molecule
    from pyPept.peptide import Peptide, serialize
    from pyPept.sequence import Sequence
    from pyPept.structure import compare_structures

    sequence, _ = parse_source(source, track_source=True)
    peptide = Peptide.from_sequence(sequence)
    emission = serialize(peptide, notation=notation)
    converted = Sequence(emission.text, warning_sink=lambda message: None)
    remapped = Peptide.from_sequence(converted, emission.occurrence_order)
    if set(peptide.connections) != set(remapped.connections):
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
    original = Molecule(sequence, depiction=None).get_molecule(fmt="ROMol")
    product = Molecule(converted, depiction=None).get_molecule(fmt="ROMol")
    if not compare_structures(original, product).exact:
        raise ValueError("Notation conversion would change the molecular structure")
    return emission.text
