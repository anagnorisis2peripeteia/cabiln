#!/usr/bin/env python3
"""Scan monomer library for structural duplicates (same canonical SMILES post-LG removal)."""
import sys, os
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers loaded\n")


def get_core_smiles(mol):
    """Get canonical SMILES with dummy atoms removed (core structure only)."""
    if mol is None:
        return None
    emol = Chem.RWMol(mol)
    # Find and remove dummy atoms (atomic num 0)
    dummies = [a.GetIdx() for a in emol.GetAtoms() if a.GetAtomicNum() == 0]
    for idx in sorted(dummies, reverse=True):
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except:
        return None


def get_full_smiles(mol):
    """Get canonical SMILES including dummies."""
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


# Group by core SMILES (post-LG, post-dummy removal)
core_groups = defaultdict(list)
full_groups = defaultdict(list)

for sym in df.index:
    mol = df.loc[sym, 'ROMol']
    core = get_core_smiles(mol)
    full = get_full_smiles(mol)
    if core:
        core_groups[core].append(sym)
    if full:
        full_groups[full].append(sym)

# Report exact duplicates (same full SMILES including dummies)
print("=" * 80)
print("EXACT DUPLICATES (identical SMILES including R-group positions)")
print("=" * 80)
exact_dupes = 0
for smi, syms in sorted(full_groups.items(), key=lambda x: -len(x[1])):
    if len(syms) > 1:
        exact_dupes += 1
        rgroups = []
        for s in syms:
            rg = df.loc[s, 'm_Rgroups'] if 'm_Rgroups' in df.columns else '?'
            mtype = df.loc[s, 'type'] if 'type' in df.columns else '?'
            rgroups.append(f"    {s:20s} type={mtype}  LG={rg}")
        print(f"\n  SMILES: {smi}")
        print("\n".join(rgroups))

print(f"\n  Total exact duplicate groups: {exact_dupes}")

# Report core duplicates (same structure after removing dummies)
print("\n\n" + "=" * 80)
print("CORE DUPLICATES (same structure after R-group/dummy removal)")
print("These may differ only in which atom the dummy is attached to")
print("=" * 80)
core_dupes = 0
for core_smi, syms in sorted(core_groups.items(), key=lambda x: -len(x[1])):
    if len(syms) > 1:
        # Skip if already reported as exact duplicate
        full_smis = set()
        for s in syms:
            mol = df.loc[s, 'ROMol']
            full_smis.add(get_full_smiles(mol))
        if len(full_smis) == 1:
            continue  # already covered by exact dupes

        core_dupes += 1
        print(f"\n  Core SMILES: {core_smi}")
        for s in syms:
            mol = df.loc[s, 'ROMol']
            full = get_full_smiles(mol)
            rg = df.loc[s, 'm_Rgroups'] if 'm_Rgroups' in df.columns else '?'
            mtype = df.loc[s, 'type'] if 'type' in df.columns else '?'
            print(f"    {s:20s} type={mtype:6s} full={full}")
            print(f"    {'':20s} LG={rg}")

print(f"\n  Total core-only duplicate groups: {core_dupes}")

# Case-insensitive name collisions
print("\n\n" + "=" * 80)
print("CASE-INSENSITIVE NAME COLLISIONS")
print("=" * 80)
name_groups = defaultdict(list)
for sym in df.index:
    name_groups[sym.lower()].append(sym)

case_dupes = 0
for lower, syms in sorted(name_groups.items()):
    if len(syms) > 1:
        case_dupes += 1
        print(f"\n  {lower}:")
        for s in syms:
            mol = df.loc[s, 'ROMol']
            full = get_full_smiles(mol)
            rg = df.loc[s, 'm_Rgroups'] if 'm_Rgroups' in df.columns else '?'
            mtype = df.loc[s, 'type'] if 'type' in df.columns else '?'
            print(f"    {s:20s} type={mtype:6s} SMILES={full}")
            print(f"    {'':20s} LG={rg}")

print(f"\n  Total case collisions: {case_dupes}")

# Summary
print("\n\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"  Total monomers: {len(df)}")
print(f"  Exact duplicate groups: {exact_dupes}")
print(f"  Core-only duplicate groups: {core_dupes}")
print(f"  Case-insensitive name collisions: {case_dupes}")
