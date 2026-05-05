#!/usr/bin/env python3
"""
Phase 3: Update SDF m_Rgroups and mol objects to match what pre_activate
produces from the fully-restored SMILES. This ensures round-trip consistency.

Skips the 6 hydrocarbon caps that pre_activate can't handle.
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers\n")

SKIP_CAPS = {'Dmb', 'Ph', '_Bn', '_Et', '_Me', 'otbu'}

LG_ATOMS = {
    '[H]': (1, 0), '[OH]': (8, 1),
    '[Cl]': (17, 0), '[Br]': (35, 0), '[I]': (53, 0),
}

def restore_full_smiles(mol, rgroups):
    emol = Chem.RWMol(Chem.RWMol(mol))
    try:
        Chem.Kekulize(emol, clearAromaticFlags=True)
    except Exception:
        pass
    dummies = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummies.append((atom.GetIdx(), atom.GetIsotope()))
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
            replacements.append((dummy_idx, 'replace_halide', LG_ATOMS[lg][0]))
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
    except Exception:
        return None

def parse_rgroups(rg):
    if isinstance(rg, list):
        return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]

def canonical(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else None

def format_rgroups(leaving_dict, max_slot=6):
    parts = []
    for i in range(1, max_slot + 1):
        if i in leaving_dict:
            parts.append(leaving_dict[i])
        else:
            parts.append('None')
    while parts and parts[-1] == 'None':
        pass  # keep trailing Nones for fixed-width format
    return ','.join(parts)

updated = []
skipped_caps = []
failed_restore = []
failed_activate = []
already_match = []

for sym in sorted(df.index):
    if sym in SKIP_CAPS:
        skipped_caps.append(sym)
        continue

    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue

    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    stored_chuckles = Chem.MolToSmiles(mol)

    full_smi = restore_full_smiles(mol, rgroups)
    if full_smi is None:
        failed_restore.append(sym)
        continue

    try:
        result = pre_activate(full_smi)
    except Exception as e:
        failed_activate.append((sym, str(e)[:60]))
        continue

    stored_canon = canonical(stored_chuckles)
    result_canon = canonical(result.chuckles)

    stored_lgs = {}
    for i, lg in enumerate(rgroups, 1):
        if lg is not None:
            stored_lgs[i] = lg

    if stored_canon == result_canon and result.leaving == stored_lgs:
        already_match.append(sym)
        continue

    # Mismatch — update SDF entry
    new_mol = Chem.MolFromSmiles(result.chuckles)
    if new_mol is None:
        print(f"  WARNING: cannot parse result CHUCKLES for {sym}: {result.chuckles}")
        continue

    AllChem.Compute2DCoords(new_mol)

    # Build m_Rgroups string
    max_slot = max(result.leaving.keys()) if result.leaving else 0
    max_slot = max(max_slot, 6)  # keep minimum 6 for compatibility
    rg_parts = []
    for i in range(1, max_slot + 1):
        if i in result.leaving:
            rg_parts.append(result.leaving[i])
        else:
            rg_parts.append('None')
    new_rgroups = ','.join(rg_parts)

    df.at[sym, 'ROMol'] = new_mol
    df.at[sym, 'm_Rgroups'] = new_rgroups
    updated.append(sym)

print(f"Already match:    {len(already_match)}")
print(f"Updated:          {len(updated)}")
print(f"Skipped caps:     {len(skipped_caps)}")
print(f"Failed restore:   {len(failed_restore)}")
print(f"Failed activate:  {len(failed_activate)}")
print(f"Total:            {len(df)}")

if updated:
    print(f"\n--- Updated monomers ({len(updated)}) ---")
    for sym in updated:
        rg = df.loc[sym, 'm_Rgroups']
        smi = Chem.MolToSmiles(df.loc[sym, 'ROMol'])
        print(f"  {sym:20s} {smi[:50]}  rg={rg[:40]}")

if failed_activate:
    print(f"\n--- Failed activate ({len(failed_activate)}) ---")
    for sym, err in failed_activate:
        print(f"  {sym}: {err}")

# Save
print(f"\nSaving SDF...")
df_save = df.copy().reset_index()
df_save = df_save.rename(columns={'index': 'symbol'})
PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print("  Done")
