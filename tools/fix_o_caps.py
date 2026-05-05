#!/usr/bin/env python3
"""
Fix monomers with [2*]OC(=O) convention that breaks restore.
These have dummy on O with LG=[OH], which restores to peroxide OO.
Change to [2*]C(=O)O (dummy on C) or fix LG to [H].

Also fixes any remaining [2*]OC(=O) diacid/fatty acid caps.
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')

def parse_rgroups(rg):
    if isinstance(rg, list):
        return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]

# Find all monomers with [2*]O pattern (dummy on oxygen of ester/acid)
broken = []
for sym in sorted(df.index):
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue
    smi = Chem.MolToSmiles(mol)
    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])

    # Check for [n*]O pattern where dummy bonds to O (not O bonding to C of COOH)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso < 1 or iso > len(rgroups):
            continue
        lg = rgroups[iso - 1]
        if lg != '[OH]':
            continue
        # Dummy with LG=[OH] — check if it bonds to an O atom
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                broken.append((sym, smi, iso, lg))
                break

print(f"Found {len(broken)} monomers with dummy-O bond + LG=[OH] (peroxide on restore)")
for sym, smi, slot, lg in broken:
    print(f"  {sym:20s} slot={slot} LG={lg}  {smi[:60]}")

# Fix strategy: these are all C-terminal caps where [2*] bonds to O of -OC(=O)-.
# The correct convention: [2*] should be on C of C(=O), LG=[OH].
# But we can't just swap atoms in the existing mol.
# Instead, restore the full SMILES correctly and re-process.

# For these monomers, the full SMILES is known:
# [2*]OC(=O)R with LG=[OH] → full SMILES should be HOOC-R (carboxylic acid)
# Not HO-OC(=O)R (peroxide!)

# So the LG is wrong. The correct LG for dummy-on-O is [H] (just H on the oxygen).
# OR we change the dummy position to be on C.

# Since the pipeline puts dummy on C, let's fix these by changing the CHUCKLES
# to put dummy on C with LG=[OH] (matching pipeline convention).

fixed = 0
for sym, smi, slot, lg in broken:
    mol = df.loc[sym, 'ROMol']
    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])

    # Rebuild the mol with dummy on the C instead of the O
    emol = Chem.RWMol(Chem.RWMol(mol))

    # Find the dummy atom for this slot
    dummy_idx = None
    o_idx = None
    c_idx = None
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            dummy_idx = atom.GetIdx()
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    o_idx = nb.GetIdx()
                    # Find the C bonded to O (the carbonyl carbon)
                    for nb2 in nb.GetNeighbors():
                        if nb2.GetAtomicNum() == 6 and nb2.GetIdx() != dummy_idx:
                            c_idx = nb2.GetIdx()
                    break
            break

    if dummy_idx is None or o_idx is None or c_idx is None:
        print(f"  WARNING: couldn't find pattern for {sym}, skipping")
        continue

    # Remove bond dummy-O, add bond dummy-C, remove bond O-C, add bond O stays as part of COOH
    # Actually, simpler: remove the dummy, put a new dummy on the C
    # But manipulating bonds is error-prone. Let's just generate new SMILES.

    # Strategy: replace [2*]OC(=O) with [2*]C(=O)O in the SMILES
    # and keep LG=[OH]

    # More robust: parse the SMILES, do the swap
    new_smi = smi
    # The pattern is: [n*]OC where n is the slot number
    pat = f'[{slot}*]O'
    if pat not in new_smi:
        # Try canonical form
        pat = f'[{slot}*]OC(=O)'
        replacement = f'[{slot}*]C(=O)O'
        if pat in new_smi:
            new_smi = new_smi.replace(pat, replacement)
        else:
            print(f"  WARNING: pattern not found for {sym}: {smi}")
            continue
    else:
        # Replace [n*]O with the O being moved to after C(=O)
        # This is complex due to SMILES ordering. Let's use the full SMILES approach.
        # Restore fully first, then note we need to fix
        pass

    if new_smi == smi:
        # Pattern replacement didn't work, use the mol manipulation approach
        # Remove dummy-O bond, remove O-C bond, add dummy-C bond, add O back
        try:
            emol.RemoveBond(dummy_idx, o_idx)
            emol.RemoveBond(o_idx, c_idx)
            emol.AddBond(dummy_idx, c_idx, Chem.BondType.SINGLE)
            emol.AddBond(o_idx, c_idx, Chem.BondType.SINGLE)
            # O now has one less neighbor (was: dummy+C, now: C only)
            # Set explicit Hs
            emol.GetAtomWithIdx(o_idx).SetNumExplicitHs(1)
            Chem.SanitizeMol(emol)
            new_smi = Chem.MolToSmiles(emol)
        except Exception as e:
            print(f"  WARNING: mol manipulation failed for {sym}: {e}")
            continue

    new_mol = Chem.MolFromSmiles(new_smi)
    if new_mol is None:
        print(f"  WARNING: couldn't parse new SMILES for {sym}: {new_smi}")
        continue

    AllChem.Compute2DCoords(new_mol)
    df.at[sym, 'ROMol'] = new_mol
    print(f"  FIXED {sym}: {smi[:40]} -> {Chem.MolToSmiles(new_mol)[:40]}")
    fixed += 1

print(f"\nFixed {fixed}/{len(broken)} monomers")

if fixed > 0:
    print(f"\nSaving SDF...")
    df_save = df.copy().reset_index()
    df_save = df_save.rename(columns={'index': 'symbol'})
    PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                         properties=list(df_save.columns.drop(['ROMol'])))
    print("  Done")
