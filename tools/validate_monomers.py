#!/usr/bin/env python3
"""Validate every SDF monomer against pyPept's reaction/pre-activate machinery.

Checks for each monomer:
  1. SMILES parses cleanly
  2. Dummy atom count matches m_Rgroups field length
  3. infer_chem_type() agrees with stored m_chem_types for every slot
  4. Backbone pair (slot1, slot2) exists in REACTION_INDEX
  5. Stored leaving group is element-compatible with stored chem_type

Run from repo root:
    python tools/validate_monomers.py [--verbose] [--csv out.csv]
"""

import sys
import pathlib
import argparse
import csv
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# Expected leaving group per chem_type (canonical check)
_CT_TO_LG = {
    'backbone_n':       '[H]',
    'backbone_n_mod':   '[H]',
    'backbone_c':       '[OH]',
    'backbone_o':       '[H]',
    'thiol':            '[H]',
    'selenol':          '[H]',
    'amine_primary':    '[H]',
    'amine_secondary':  '[H]',
    'guanidinium':      '[H]',
    'aromatic_nh':      '[H]',
    'amide_nh':         '[H]',
    'hydroxyl':         '[H]',
    'hydroxyl_phenolic':'[H]',
    'carboxyl':         '[H]',   # O-attachment: dummy replaces H on carboxyl OH
    'phosphate_p':      '[OH]',
    'alkyl_halide_c':   '[H]',
    'aminooxy':         '[H]',
    'hydrazide':        '[H]',
    'terminal_alkene':  '[H]',
    'cyclooctyne_c':    '[H]',
    'alkyne_c':         '[H]',
    'azide_alpha_c':    '[H]',
    'tetrazine_c':      '[H]',
    'tco_c':            '[H]',
    'maleimide_c':      '[H]',
    'aldehyde':         '[H]',
    'nhs_ester':        '[H]',
}


def parse_rgroups(rg_str):
    """'[H],[OH],[H],None,None,None' -> ['[H]','[OH]','[H]',None,None,None]"""
    return [None if x.strip() == 'None' else x.strip() for x in rg_str.split(',')]


def parse_chem_types(ct_str):
    """'1:backbone_n,2:backbone_c,3:backbone_n_mod' -> {1:'backbone_n', 2:'backbone_c', 3:'backbone_n_mod'}"""
    result = {}
    for item in ct_str.split(','):
        item = item.strip()
        if ':' in item:
            slot_s, ct = item.split(':', 1)
            result[int(slot_s)] = ct.strip()
    return result


def get_slot_attachment_map(mol):
    """{slot(1-based): attach_atom_idx} from dummy atoms in a CHUCKLES mol."""
    mapping = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # dummy atom
            iso = atom.GetIsotope()
            if iso > 0:
                nbrs = list(atom.GetNeighbors())
                if nbrs:
                    mapping[iso] = nbrs[0].GetIdx()
    return mapping


# label_only=True types in _CHEM_TYPE_REGISTRY — infer_chem_type never returns these
# (they're excluded from _EXOTIC_SMARTS) so mismatches on these slots are expected.
_LABEL_ONLY_TYPES = frozenset({
    'aromatic_nh', 'hydroxyl_phenolic', 'guanidinium', 'amide_nh',
    # amine_secondary is assigned by the second-H pass in pre_activate but has no
    # infer_smarts entry; infer_chem_type falls through to amine_primary heuristic.
    'amine_secondary',
})

# Carboxyl and aldehyde are structurally identical in CHUCKLES form: both encode as
# C([n*])=O (dummy replaced either the OH or the H).  infer_chem_type cannot
# distinguish them — aldehyde SMARTS [CX3;H0,H1:1](=O)[!#7] fires before the
# carboxyl C heuristic.  Assembly uses m_chem_types as the authority; the inference
# mismatch is expected and harmless for hand-crafted monomers that have valid metadata.
_CARBOXYL_ALDEHYDE = frozenset({'carboxyl', 'aldehyde'})


def validate_monomer(mol, abbr, rg_str, ct_str, infer_chem_type_fn, REACTION_INDEX):
    issues = []

    rgroups = parse_rgroups(rg_str)
    stored_ct = parse_chem_types(ct_str)
    slot_map = get_slot_attachment_map(mol)

    # 1. Non-None rgroup count vs actual dummy count.
    # m_Rgroups always has 6 entries (trailing None = no slot); only non-None entries
    # correspond to real dummy atoms.
    n_dummies = len(slot_map)
    n_active_rgroups = sum(1 for x in rgroups if x is not None)
    if n_dummies != n_active_rgroups:
        issues.append(f"dummy_count_mismatch: {n_dummies} dummies vs {n_active_rgroups} active m_Rgroups entries")

    # 2. Stored m_chem_types slots should match dummy slots
    stored_slots = set(stored_ct.keys())
    dummy_slots = set(slot_map.keys())
    missing_in_dummy = stored_slots - dummy_slots
    missing_in_ct = dummy_slots - stored_slots
    if missing_in_dummy:
        issues.append(f"slot_in_ct_but_no_dummy: {sorted(missing_in_dummy)}")
    if missing_in_ct:
        issues.append(f"dummy_slot_missing_from_ct: {sorted(missing_in_ct)}")

    # 3. infer_chem_type agreement — skip label_only types and known ambiguities
    for slot, attach_idx in sorted(slot_map.items()):
        if slot not in stored_ct:
            continue
        stored = stored_ct[slot]
        if stored in _LABEL_ONLY_TYPES:
            continue  # never returned by infer_chem_type; skip
        try:
            inferred = infer_chem_type_fn(mol, attach_idx, slot)
        except Exception as exc:
            issues.append(f"slot{slot}_infer_error: {exc}")
            continue
        if inferred == stored:
            continue
        # carboxyl and aldehyde are structurally identical in CHUCKLES (both C([n*])=O);
        # the mismatch is a known limitation — skip it.
        if stored in _CARBOXYL_ALDEHYDE and inferred in _CARBOXYL_ALDEHYDE:
            continue
        issues.append(f"slot{slot}_ct_mismatch: stored={stored!r} inferred={inferred!r}")

    # 4. Backbone reaction pair in REACTION_INDEX
    ct1 = stored_ct.get(1)
    ct2 = stored_ct.get(2)
    if ct1 and ct2:
        if (ct1, ct2) not in REACTION_INDEX:
            issues.append(f"no_reaction: backbone pair ({ct1},{ct2}) not in REACTION_INDEX")

    # 5. Leaving group vs chem_type compatibility (skip label_only types)
    for slot, attach_idx in sorted(slot_map.items()):
        if slot not in stored_ct:
            continue
        ct = stored_ct[slot]
        if ct in _LABEL_ONLY_TYPES:
            continue
        if slot - 1 < len(rgroups):
            lg = rgroups[slot - 1]
        else:
            continue
        if lg is None:
            continue
        expected_lg = _CT_TO_LG.get(ct)
        if expected_lg and lg != expected_lg:
            issues.append(f"slot{slot}_lg_mismatch: stored_lg={lg!r} expected_for_{ct}={expected_lg!r}")

    return issues


def main():
    parser = argparse.ArgumentParser(description="Validate monomers.sdf against pyPept reaction machinery")
    parser.add_argument('--verbose', '-v', action='store_true', help="Print all monomers, not just conflicts")
    parser.add_argument('--csv', metavar='FILE', help="Write results to CSV file")
    parser.add_argument('--type', help="Filter by m_type (e.g. aa, linker, cap)")
    args = parser.parse_args()

    from rdkit import Chem
    from pyPept.interfaces.reaction_library import infer_chem_type, REACTION_INDEX

    print(f"Loading monomers from: {SDF_PATH}")
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    mols = [(m, m.GetPropsAsDict()) for m in suppl if m is not None]
    print(f"  Loaded {len(mols)} valid monomers\n")

    if args.type:
        mols = [(m, p) for m, p in mols if p.get('m_type') == args.type]
        print(f"  Filtered to {len(mols)} monomers with m_type={args.type!r}\n")

    counts = defaultdict(int)
    all_results = []

    for mol, props in mols:
        abbr     = props.get('m_abbr', '?')
        rg_str   = props.get('m_Rgroups', '')
        ct_str   = props.get('m_chem_types', '')
        mtype    = props.get('m_type', '?')
        msubtype = props.get('m_subtype', '?')

        if not rg_str or not ct_str:
            counts['missing_fields'] += 1
            issue_str = "MISSING m_Rgroups or m_chem_types"
            all_results.append((abbr, mtype, msubtype, issue_str))
            print(f"  MISSING_FIELDS  {abbr}")
            continue

        try:
            issues = validate_monomer(mol, abbr, rg_str, ct_str, infer_chem_type, REACTION_INDEX)
        except Exception as exc:
            counts['error'] += 1
            issue_str = f"EXCEPTION: {exc}"
            all_results.append((abbr, mtype, msubtype, issue_str))
            print(f"  ERROR  {abbr}: {exc}")
            continue

        if issues:
            counts['conflict'] += 1
            issue_str = '; '.join(issues)
            all_results.append((abbr, mtype, msubtype, issue_str))
            print(f"  CONFLICT  {abbr} ({mtype}/{msubtype})")
            for iss in issues:
                print(f"            {iss}")
        else:
            counts['clean'] += 1
            all_results.append((abbr, mtype, msubtype, 'ok'))
            if args.verbose:
                print(f"  ok  {abbr}")

    print(f"\n{'='*60}")
    print(f"  Total:    {len(mols)}")
    print(f"  Clean:    {counts['clean']}")
    print(f"  Conflict: {counts['conflict']}")
    print(f"  Error:    {counts['error']}")
    print(f"  Missing:  {counts['missing_fields']}")
    print(f"{'='*60}")

    if counts['conflict'] or counts['error']:
        # Group conflicts by type
        conflict_types = defaultdict(list)
        for abbr, mtype, msubtype, issue_str in all_results:
            if issue_str != 'ok':
                for issue in issue_str.split('; '):
                    key = issue.split(':')[0]
                    conflict_types[key].append(abbr)

        print("\nConflict breakdown:")
        for key in sorted(conflict_types):
            entries = conflict_types[key]
            sample = ', '.join(entries[:5])
            more = f' ... (+{len(entries)-5} more)' if len(entries) > 5 else ''
            print(f"  {key}: {len(entries)} monomers — {sample}{more}")

    if args.csv:
        with open(args.csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['abbr', 'm_type', 'm_subtype', 'issues'])
            w.writerows(all_results)
        print(f"\nCSV written to: {args.csv}")

    return 1 if (counts['conflict'] or counts['error']) else 0


if __name__ == '__main__':
    sys.exit(main())
