"""Restore unbound attachment dummies for assembly and monomer depictions."""

from __future__ import annotations

from collections.abc import Mapping
import re

from rdkit import Chem


def leaving_family(value):
    """Isotopic variants retain the same attachment chemistry."""
    return re.sub(r'^\[\d+(H|OH|Cl|Br|I)\]$', r'[\1]', value) if value else value


def restore_leaving_groups(
    mol: Chem.Mol,
    leaving_by_isotope: Mapping[int, str | None],
    *,
    sanitize: bool = True,
    slots=None,
) -> Chem.Mol:
    """Return a copy with each dummy replaced by its single-atom leaving group.

    Neutral, unlabeled hydrogen is implicit; isotope-labeled hydrogen remains
    explicit. Charge, isotopes and residue metadata survive replacement. As in
    the existing assembly path, missing leaving metadata means implicit H.
    Multi-atom leaving groups are unsupported and fail instead of truncating.
    Assembly may defer sanitization to its existing contextual error handler.
    ``slots`` optionally restores only selected dummies, retaining the others.
    """
    editable = Chem.RWMol(mol)
    try:
        Chem.Kekulize(editable, clearAromaticFlags=True)
    except ValueError:
        # Some intermediates only become sanitizable after dummy removal.
        # The final sanitization still reports an invalid restored structure.
        pass

    to_remove = []
    for atom in editable.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        if slots is not None and atom.GetIsotope() not in slots:
            continue
        leaving = leaving_by_isotope.get(atom.GetIsotope())
        if isinstance(leaving, str):
            leaving = leaving.strip()
            if leaving.lower() in ("", "none"):
                leaving = None
        source = None
        if leaving is not None:
            parser = Chem.SmilesParserParams()
            parser.removeHs = False
            parsed = Chem.MolFromSmiles(leaving, parser)
            if (
                parsed is None
                or parsed.GetNumAtoms() != 1
                or parsed.GetAtomWithIdx(0).GetAtomicNum() == 0
            ):
                raise ValueError(
                    f"Unsupported leaving group {leaving!r}: terminal restoration "
                    "requires a single atom."
                )
            source = parsed.GetAtomWithIdx(0)

        ordinary_hydrogen = source is None or (
            source.GetAtomicNum() == 1
            and source.GetIsotope() == 0
            and source.GetFormalCharge() == 0
        )
        if ordinary_hydrogen:
            for neighbor in atom.GetNeighbors():
                neighbor.SetNoImplicit(False)
            to_remove.append(atom.GetIdx())
            continue

        # Recalculate valence in the bonded product. Copying the free [OH]
        # fragment's radical state, for example, would create a radical acid.
        replacement = Chem.Atom(source.GetAtomicNum())
        replacement.SetFormalCharge(source.GetFormalCharge())
        replacement.SetIsotope(source.GetIsotope())
        if atom.GetPDBResidueInfo() is not None:
            replacement.SetMonomerInfo(atom.GetPDBResidueInfo())
        editable.ReplaceAtom(atom.GetIdx(), replacement, preserveProps=True)

    for index in reversed(to_remove):
        editable.RemoveAtom(index)
    restored = editable.GetMol()
    if sanitize:
        try:
            Chem.SanitizeMol(restored)
        except ValueError as exc:
            raise ValueError(
                f"Leaving-group restoration produced an invalid structure: {exc}"
            ) from exc
    return restored
