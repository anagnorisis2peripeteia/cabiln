#!/usr/bin/env python3
"""Check and fix all fatty acid / diacid chain lengths."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')

# All fatty acid / diacid caps
fa_caps = ['C14FA', 'C16FA', 'C18FA', 'C20FA', 'C22FA',
           'C14dAcid', 'C16dAcid', 'C22dAcid', 'Dod_dAcid',
           'Succ', 'Glt', 'Adp', 'Pim', 'Sub', 'Azel', 'Seb',
           'Myr', 'Pal', 'Ste', 'Ole']

print("\nCurrent chain lengths:")
print(f"{'Symbol':15s} {'Carbons':8s} {'Expected':10s} {'SMILES'}")
print("-" * 90)

for sym in fa_caps:
    if sym not in df.index:
        print(f"  {sym}: NOT IN SDF")
        continue
    mol = df.loc[sym, 'ROMol']
    smi = Chem.MolToSmiles(mol) if mol else '?'
    # Count carbons in the SMILES (excluding the dummy)
    if mol:
        n_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    else:
        n_carbons = 0
    print(f"  {sym:15s} {n_carbons:3d}C     {smi[:65]}")

# Expected carbon counts for diacid caps:
# CnFA/CndAcid = n-carbon diacid, one COOH as R2, other as free COOH
# After R2 dummy, SMILES = [2*]C(=O)(CH2)_{n-2}C(=O)O
# Total C in molecule = n (the diacid has n carbons)
# Dod_dAcid = C12 (dodecanedioic)
# Succ = C4, Glt = C5, Adp = C6, Pim = C7, Sub = C8, Azel = C9, Seb = C10
# Mono-acids: Myr = C14, Pal = C16, Ste = C18, Ole = C18 (unsaturated)

EXPECTED = {
    'C14FA': 14, 'C16FA': 16, 'C18FA': 18, 'C20FA': 20, 'C22FA': 22,
    'C14dAcid': 14, 'C16dAcid': 16, 'C22dAcid': 22,
    'Dod_dAcid': 12,
    'Succ': 4, 'Glt': 5, 'Adp': 6, 'Pim': 7, 'Sub': 8, 'Azel': 9, 'Seb': 10,
    'Myr': 14, 'Pal': 16, 'Ste': 18, 'Ole': 18,
}

print("\n\nVERIFICATION:")
print(f"{'Symbol':15s} {'Actual':8s} {'Expected':10s} {'Status'}")
print("-" * 50)

fixes_needed = {}
for sym, expected_c in sorted(EXPECTED.items()):
    if sym not in df.index:
        print(f"  {sym:15s} {'N/A':8s} {expected_c:3d}C       NOT IN SDF")
        continue
    mol = df.loc[sym, 'ROMol']
    actual_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    status = "OK" if actual_c == expected_c else f"WRONG (off by {actual_c - expected_c})"
    print(f"  {sym:15s} {actual_c:3d}C     {expected_c:3d}C       {status}")
    if actual_c != expected_c:
        fixes_needed[sym] = expected_c

if fixes_needed:
    print(f"\n\nFIXING {len(fixes_needed)} entries...")
    for sym, expected_c in fixes_needed.items():
        mol = df.loc[sym, 'ROMol']
        smi = Chem.MolToSmiles(mol)
        # Determine if this is a diacid (has free COOH end) or mono-acid
        is_diacid = 'C(=O)O' in smi or '(=O)O' in smi.replace('[2*]C(=O)', '', 1)

        if is_diacid:
            # Diacid: [2*]C(=O) + (n-2) CH2 + C(=O)O
            n_ch2 = expected_c - 2
            new_smi = '[2*]C(=O)' + 'C' * n_ch2 + 'C(=O)O'
        else:
            # Mono-acid (fatty acid): [2*]C(=O) + (n-1) CH2 + terminal CH3
            # Actually mono fatty acid caps: [2*]C(=O)(CH2)_{n-2}CH3
            # n carbons total in the acyl chain including carbonyl C
            n_ch2 = expected_c - 2  # -1 for carbonyl C, -1 for terminal CH3
            new_smi = '[2*]C(=O)' + 'C' * n_ch2 + 'C'

        new_mol = Chem.MolFromSmiles(new_smi)
        if new_mol is None:
            print(f"  ERROR: can't parse {new_smi} for {sym}")
            continue

        # Verify carbon count
        new_c = sum(1 for a in new_mol.GetAtoms() if a.GetAtomicNum() == 6)
        if new_c != expected_c:
            print(f"  ERROR: generated SMILES has {new_c}C, expected {expected_c}C for {sym}")
            continue

        AllChem.Compute2DCoords(new_mol)
        df.at[sym, 'ROMol'] = new_mol
        print(f"  FIXED {sym}: {smi[:50]} -> {Chem.MolToSmiles(new_mol)[:50]}")

    # Save
    print("\nSaving SDF...")
    df_save = df.copy().reset_index()
    df_save = df_save.rename(columns={'index': 'symbol'})
    PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                         properties=list(df_save.columns.drop(['ROMol'])))
    print("  Done")
else:
    print("\n  All chain lengths correct!")
