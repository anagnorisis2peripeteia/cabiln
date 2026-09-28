"""Inspect the attachment chemistry used by assembly and the dynamic palette."""

from __future__ import annotations

from typing import NamedTuple

from pyPept.interfaces.reaction_library import REACTION_INDEX, infer_chem_type


class ResolvedConnection(NamedTuple):
    chem_type1: str
    chem_type2: str
    reaction: dict | None


def reaction_for_types(chem_type1, chem_type2):
    """Return the supported reaction; the index already includes both orders."""
    return REACTION_INDEX.get((chem_type1, chem_type2))


def _site_metadata(mol, leaving_groups):
    if leaving_groups is None:
        leaving_groups = mol.GetProp("m_Rgroups") if mol.HasProp("m_Rgroups") else ""
    if isinstance(leaving_groups, str):
        leaving_groups = leaving_groups.split(",")
    declared = {}
    if mol.HasProp("m_chem_types"):
        for field in mol.GetProp("m_chem_types").split(","):
            slot, separator, name = field.partition(":")
            if separator and slot.strip().isdigit():
                declared[int(slot)] = name.strip()
    return leaving_groups, declared


def _describe_site(mol, atom, leaving_groups, declared):
    slot = atom.GetIsotope()
    if slot < 1 or atom.GetDegree() != 1:
        raise ValueError("Attachment dummies need a positive slot and one neighbour")
    leaving = leaving_groups[slot - 1] if slot <= len(leaving_groups) else None
    leaving = leaving.strip() if isinstance(leaving, str) else leaving
    if leaving in ("", "None", "none"):
        leaving = None
    return {
        "slot": slot,
        "chem_type": infer_chem_type(
            mol, atom.GetNeighbors()[0].GetIdx(), slot=slot, leaving=leaving
        ),
        "leaving": leaving or "",
        "declared_chem_type": declared.get(slot, ""),
    }


def attachment_site(mol, slot, leaving_groups=None):
    """Resolve only the selected numbered site, without enriching other sites."""
    atoms = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot
    ]
    if not atoms:
        raise ValueError(f"Monomer has no R{slot} attachment point")
    if len(atoms) != 1:
        raise ValueError(f"Monomer has multiple R{slot} attachment points")
    return _describe_site(mol, atoms[0], *_site_metadata(mol, leaving_groups))


def resolve_connection(
    mol1, slot1, mol2, slot2, *, leaving_groups1=None, leaving_groups2=None
):
    """Resolve effective site types and their supported reaction in endpoint order.

    A missing reaction is a chemistry incompatibility, distinct from a missing
    attachment site. Callers add their occurrence and source context to errors.
    Reaction execution retains responsibility for matching template orientation.
    """
    type1 = attachment_site(mol1, slot1, leaving_groups1)["chem_type"]
    type2 = attachment_site(mol2, slot2, leaving_groups2)["chem_type"]
    return ResolvedConnection(type1, type2, reaction_for_types(type1, type2))


def attachment_sites(mol, leaving_groups=None):
    """Describe numbered dummy sites using the same detection as assembly.

    Declared labels remain available as metadata. Reaction compatibility uses
    the effective chemistry type inferred from the structure, slot and leaving
    group, including the assembly engine's supported declaration overrides.
    """
    leaving_groups, declared = _site_metadata(mol, leaving_groups)
    sites = [
        _describe_site(mol, atom, leaving_groups, declared)
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0
    ]
    return sorted(sites, key=lambda site: site["slot"])
