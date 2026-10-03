"""Inspect the attachment chemistry used by assembly and the dynamic palette."""

from __future__ import annotations

from typing import NamedTuple

from pyPept.interfaces.reaction_library import REACTION_INDEX
from pyPept.site_chemistry import Perception

_SUPPORTED_TYPES = frozenset(kind for pair in REACTION_INDEX for kind in pair)


def _attachment_idx(mol, slot):
    """Return attachment atom index for R-group slot (1-based), or None.

    Locates the dummy atom whose isotope equals slot (CHUCKLES convention:
    slot 1 → [1*], slot 2 → [2*], …) and returns its neighbour's index —
    i.e. the heavy atom where the inter-monomer bond will form.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            nb = atom.GetNeighbors()
            return nb[0].GetIdx() if nb else None
    return None


def _rgroup_atom_idx(mol, slot):
    """Return the dummy atom index for R-group slot (1-based), or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            return atom.GetIdx()
    return None


def _slot_for_attachment(mol, atom_idx):
    """Return R-group slot (1-based) whose dummy neighbours atom_idx, or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            for nb in atom.GetNeighbors():
                if nb.GetIdx() == atom_idx:
                    return atom.GetIsotope()
    return None


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


def _describe_site(perception, atom, leaving_groups, declared):
    slot = atom.GetIsotope()
    if slot < 1 or atom.GetDegree() != 1:
        raise ValueError("Attachment dummies need a positive slot and one neighbour")
    leaving = leaving_groups[slot - 1] if slot <= len(leaving_groups) else None
    leaving = leaving.strip() if isinstance(leaving, str) else leaving
    if leaving in ("", "None", "none"):
        leaving = None
    chemistry = perception.classify(
        atom.GetNeighbors()[0].GetIdx(), slot, leaving, declared.get(slot)
    ).reaction_type
    return {
        "slot": slot,
        "chem_type": chemistry,
        "supported": chemistry in _SUPPORTED_TYPES,
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
    return _describe_site(Perception(mol), atoms[0], *_site_metadata(mol, leaving_groups))


def declaration_is_compatible(mol, site):
    """Compare a declaration with resolved chemistry and documented site roles.

    Historical role labels are equivalent only in the corresponding structural
    context. In particular, an amide or guanidine never acquires the primary-
    amine second-H alias just because its old declaration says ``secondary``.
    """
    declared = site['declared_chem_type']
    effective = site['chem_type']
    if declared == effective:
        return True
    dummy = next(atom for atom in mol.GetAtoms()
                 if atom.GetAtomicNum() == 0 and atom.GetIsotope() == site['slot'])
    atom = dummy.GetNeighbors()[0]
    if (declared == 'hydroxyl_phenolic' and effective == 'aryl_phenol_o'
            and atom.GetAtomicNum() == 8):
        return any(nb.GetAtomicNum() == 6 and nb.GetIsAromatic()
                   for nb in atom.GetNeighbors())
    if (declared == 'amine_secondary' and effective == 'amine_primary'
            and atom.GetAtomicNum() == 7):
        slots = sorted(nb.GetIsotope() for nb in atom.GetNeighbors()
                       if nb.GetAtomicNum() == 0)
        return len(slots) == 2 and site['slot'] == slots[1]
    if (declared == 'backbone_o' and effective == 'hydroxyl'
            and site['slot'] == 2 and atom.GetAtomicNum() == 8):
        return any(nb.GetAtomicNum() == 6 and not nb.GetIsAromatic()
                   for nb in atom.GetNeighbors())
    return False


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
    perception = Perception(mol)
    sites = [
        _describe_site(perception, atom, leaving_groups, declared)
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0
    ]
    if len({site['slot'] for site in sites}) != len(sites):
        raise ValueError('attachment slot numbers must be unique')
    return sorted(sites, key=lambda site: site["slot"])
