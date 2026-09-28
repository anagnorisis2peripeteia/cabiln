"""Inspect the attachment chemistry used by assembly and the dynamic palette."""

from __future__ import annotations

from pyPept.interfaces.reaction_library import infer_chem_type


def attachment_sites(mol, leaving_groups=None):
    """Describe numbered dummy sites using the same detection as assembly.

    Declared labels remain available as metadata. Reaction compatibility uses
    the effective chemistry type inferred from the structure, slot and leaving
    group, including the assembly engine's supported declaration overrides.
    """
    if leaving_groups is None:
        raw = mol.GetProp("m_Rgroups") if mol.HasProp("m_Rgroups") else ""
        leaving_groups = raw.split(",")
    declared = {}
    if mol.HasProp("m_chem_types"):
        for field in mol.GetProp("m_chem_types").split(","):
            slot, separator, name = field.partition(":")
            if separator and slot.strip().isdigit():
                declared[int(slot)] = name.strip()
    sites = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        slot = atom.GetIsotope()
        if slot < 1 or atom.GetDegree() != 1:
            raise ValueError(
                "Attachment dummies need a positive slot and one neighbour"
            )
        leaving = leaving_groups[slot - 1] if slot <= len(leaving_groups) else None
        leaving = leaving.strip() if isinstance(leaving, str) else leaving
        if leaving in ("", "None", "none"):
            leaving = None
        chem_type = infer_chem_type(
            mol, atom.GetNeighbors()[0].GetIdx(), slot=slot, leaving=leaving
        )
        sites.append(
            {
                "slot": slot,
                "chem_type": chem_type,
                "leaving": leaving or "",
                "declared_chem_type": declared.get(slot, ""),
            }
        )
    return sorted(sites, key=lambda site: site["slot"])
