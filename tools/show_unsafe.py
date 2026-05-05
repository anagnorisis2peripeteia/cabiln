#!/usr/bin/env python3
"""Show the unsafe/skipped monomers from the safe update."""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')

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

def get_slot_numbers(smiles):
    return sorted(set(int(m) for m in re.findall(r'\[(\d+)\*\]', smiles)))

mismatched = []

for sym in sorted(df.index):
    if sym in SKIP_CAPS:
        continue
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue
    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    stored_chuckles = Chem.MolToSmiles(mol)
    full_smi = restore_full_smiles(mol, rgroups)
    if full_smi is None:
        continue
    try:
        result = pre_activate(full_smi)
    except Exception:
        continue
    stored_canon = canonical(stored_chuckles)
    result_canon = canonical(result.chuckles)
    stored_lgs = {}
    for i, lg in enumerate(rgroups, 1):
        if lg is not None:
            stored_lgs[i] = lg
    if stored_canon == result_canon and result.leaving == stored_lgs:
        continue
    mismatched.append((sym, stored_canon, result_canon, stored_lgs, result.leaving))

print(f"Remaining mismatches: {len(mismatched)}")
for sym, stored, result, stored_lg, result_lg in mismatched:
    s_slots = get_slot_numbers(stored)
    r_slots = get_slot_numbers(result)
    print(f"\n{sym}:")
    print(f"  stored slots: {s_slots}  result slots: {r_slots}")
    print(f"  stored: {stored[:70]}")
    print(f"  result: {result[:70]}")
    if stored_lg != result_lg:
        print(f"  stored LG: {stored_lg}")
        print(f"  result LG: {result_lg}")
