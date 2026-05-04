#!/usr/bin/env python3
"""
Full library round-trip: for every monomer in the SDF, reconstruct its full
SMILES (by restoring LGs onto the CHUCKLES), then run pre_activate and check
if the output matches the stored CHUCKLES and LGs.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from pyPept.interfaces.monomer_pipeline import pre_activate, find_cap_slots, find_backbone_slots

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers\n")

LG_ATOMS = {
    '[H]': (1, 0),    # hydrogen, valence 0 in mol context
    '[OH]': (8, 1),   # oxygen with one H
    '[Cl]': (17, 0),
    '[Br]': (35, 0),
    '[I]': (53, 0),
}

def restore_full_smiles(mol, rgroups):
    """Restore a CHUCKLES mol to full SMILES by replacing dummies with LG atoms."""
    emol = Chem.RWMol(Chem.RWMol(mol))

    # Kekulize first so aromatic NH atoms keep proper valence after dummy removal
    try:
        Chem.Kekulize(emol, clearAromaticFlags=True)
    except Exception:
        pass

    dummies = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            iso = atom.GetIsotope()
            dummies.append((atom.GetIdx(), iso))

    replacements = []
    for dummy_idx, slot in dummies:
        if slot < 1 or slot > len(rgroups):
            continue
        lg = rgroups[slot - 1]
        if lg is None:
            continue

        if lg == '[H]':
            replacements.append((dummy_idx, 'remove', None))
        elif lg == '[OH]':
            replacements.append((dummy_idx, 'replace_oh', None))
        elif lg in ('[Cl]', '[Br]', '[I]'):
            atom_num = LG_ATOMS[lg][0]
            replacements.append((dummy_idx, 'replace_halide', atom_num))
        else:
            replacements.append((dummy_idx, 'unknown', lg))

    atoms_to_remove = []
    for dummy_idx, action, param in replacements:
        if action == 'remove':
            for nb in emol.GetAtomWithIdx(dummy_idx).GetNeighbors():
                nb.SetNoImplicit(False)
            atoms_to_remove.append(dummy_idx)
        elif action == 'replace_oh':
            emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(8)
            emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
            emol.GetAtomWithIdx(dummy_idx).SetNumExplicitHs(1)
        elif action == 'replace_halide':
            emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(param)
            emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
        elif action == 'unknown':
            atoms_to_remove.append(dummy_idx)

    for idx in sorted(atoms_to_remove, reverse=True):
        emol.RemoveAtom(idx)

    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except Exception as e:
        return None


def parse_rgroups(rg):
    if isinstance(rg, list):
        return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]


def canonical(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


# Categories for results
passed = []
failed_restore = []
failed_activate = []
failed_chuckles_mismatch = []
failed_lg_mismatch = []
failed_slot_mismatch = []
skipped = []

print("=" * 80)
print("FULL LIBRARY ROUND-TRIP TEST")
print("=" * 80)

for sym in sorted(df.index):
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        skipped.append((sym, 'no mol'))
        continue

    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    stored_chuckles = Chem.MolToSmiles(mol)

    # Step 1: Restore full SMILES from CHUCKLES + LGs
    full_smi = restore_full_smiles(mol, rgroups)
    if full_smi is None:
        failed_restore.append((sym, stored_chuckles, str(rgroups)))
        continue

    # Step 2: Run pre_activate on the full SMILES
    try:
        result = pre_activate(full_smi)
    except Exception as e:
        failed_activate.append((sym, full_smi, str(e)[:80]))
        continue

    # Step 3: Compare CHUCKLES (canonical)
    stored_canon = canonical(stored_chuckles)
    result_canon = canonical(result.chuckles)

    if stored_canon != result_canon:
        failed_chuckles_mismatch.append((sym, stored_canon, result_canon, full_smi))
        continue

    # Step 4: Compare LGs
    stored_lgs = {}
    for i, lg in enumerate(rgroups, 1):
        if lg is not None:
            stored_lgs[i] = lg

    if result.leaving != stored_lgs:
        failed_lg_mismatch.append((sym, stored_lgs, result.leaving))
        continue

    passed.append(sym)

# ═══════════════════════════════════════════════════════════════
# Report
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 80}")
print("RESULTS")
print(f"{'=' * 80}")
print(f"  PASSED:              {len(passed)}")
print(f"  Failed restore:      {len(failed_restore)}")
print(f"  Failed pre_activate: {len(failed_activate)}")
print(f"  CHUCKLES mismatch:   {len(failed_chuckles_mismatch)}")
print(f"  LG mismatch:         {len(failed_lg_mismatch)}")
print(f"  Skipped:             {len(skipped)}")
print(f"  Total:               {len(df)}")

if failed_restore:
    print(f"\n--- Failed restore ({len(failed_restore)}) ---")
    for sym, chk, rg in failed_restore[:20]:
        print(f"  {sym:20s} chuckles={chk[:50]}  rg={rg[:40]}")

if failed_activate:
    print(f"\n--- Failed pre_activate ({len(failed_activate)}) ---")
    for sym, smi, err in failed_activate[:30]:
        print(f"  {sym:20s} smi={smi[:40]}  err={err}")

if failed_chuckles_mismatch:
    print(f"\n--- CHUCKLES mismatch ({len(failed_chuckles_mismatch)}) ---")
    for sym, stored, result, full_smi in failed_chuckles_mismatch[:30]:
        print(f"  {sym:20s}")
        print(f"    stored:  {stored[:60]}")
        print(f"    result:  {result[:60]}")
        print(f"    input:   {full_smi[:60]}")

if failed_lg_mismatch:
    print(f"\n--- LG mismatch ({len(failed_lg_mismatch)}) ---")
    for sym, stored, result in failed_lg_mismatch[:30]:
        print(f"  {sym:20s}")
        print(f"    stored: {stored}")
        print(f"    result: {result}")
