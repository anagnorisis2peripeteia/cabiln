"""Restore leaving groups and reagent forms for monomer depictions."""

from __future__ import annotations

from pyPept.leaving_groups import restore_leaving_groups

_CAP_REACTIONS = None


def _load_cap_reactions():
    global _CAP_REACTIONS
    if _CAP_REACTIONS is None:
        import pathlib

        import yaml

        yaml_path = (
            pathlib.Path(__file__).resolve().parents[1] / "data" / "cap_reactions.yaml"
        )
        if yaml_path.exists():
            with open(yaml_path, encoding="utf-8") as f:
                _CAP_REACTIONS = yaml.safe_load(f) or {}
        else:
            _CAP_REACTIONS = {}
    return _CAP_REACTIONS


_LG_SMILES = {
    "Cl": "[Cl]",
    "Br": "[Br]",
    "I": "[I]",
    "OH": "[OH]",
    "H": "[H]",
}


def _restore_reagent_form(mol, abbr):
    """Replace the bonding R-group dummy with the reagent's canonical leaving group.

    For an N-terminal cap (R2 bonding slot), this shows the actual reagent:
    e.g. ac with LG=Cl renders as CH₃C(=O)Cl (acetyl chloride).
    For C-terminal caps (R1 bonding slot), LG=H shows the free nucleophile.
    """
    cap_rx = _load_cap_reactions()
    info = cap_rx.get(abbr)
    if not info:
        return None, None

    reagent_lg = info.get("reagent_lg", "")
    lg_smi = _LG_SMILES.get(reagent_lg)
    if not lg_smi:
        return None, None

    props = mol.GetPropsAsDict()
    rgroups = props.get("m_Rgroups", "")
    slots = [r.strip() for r in rgroups.split(",")]

    # Determine which slot is the bonding slot for this cap type
    # R2 caps (N-terminal/electrophilic): slot 2 bonds to the amine
    # R1 caps (C-terminal/nucleophilic): slot 1 bonds to the carbonyl
    # Sidechain caps: slot 1 bonds to the sidechain
    # We detect by checking which slot has the conventional LG ([OH] or [H])
    bonding_slot = None
    if len(slots) >= 2 and slots[1] in ("[OH]", "[H]"):
        bonding_slot = 2  # R2 cap
    elif len(slots) >= 1 and slots[0] in ("[OH]", "[H]"):
        bonding_slot = 1  # R1 cap
    else:
        bonding_slot = 2  # Default to R2

    leaving_by_slot = {index + 1: leaving for index, leaving in enumerate(slots)}
    leaving_by_slot[bonding_slot] = lg_smi
    restored = restore_leaving_groups(mol, leaving_by_slot)

    reaction_type = info.get("reaction", "")
    reagent_note = info.get("reagent_note", "")
    issue = info.get("issue", "")
    meta = {
        "reaction": reaction_type,
        "reagent_lg": reagent_lg,
        "reagent_note": reagent_note,
    }
    if issue:
        meta["issue"] = issue
    return restored, meta


def _restore_leaving_groups(mol):
    """Restore a standalone monomer using the same rules as peptide assembly."""
    rgroups = mol.GetProp("m_Rgroups") if mol.HasProp("m_Rgroups") else ""
    leaving_by_slot = {
        index + 1: leaving for index, leaving in enumerate(rgroups.split(","))
    }
    return restore_leaving_groups(mol, leaving_by_slot)
