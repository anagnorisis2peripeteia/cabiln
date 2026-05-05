#!/usr/bin/env python3
"""Categorize CHUCKLES mismatches: extra slots vs slot numbering vs structural."""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')

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

def count_slots(smiles):
    return len(re.findall(r'\[\d+\*\]', smiles))

extra_slots = []
fewer_slots = []
same_slots_diff_numbering = []
same_slots_structural = []

for sym in sorted(df.index):
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
    if stored_canon == result_canon:
        continue

    stored_lgs = {}
    for i, lg in enumerate(rgroups, 1):
        if lg is not None:
            stored_lgs[i] = lg

    n_stored = count_slots(stored_canon)
    n_result = count_slots(result_canon)

    if n_result > n_stored:
        extra_slots.append((sym, n_stored, n_result, stored_canon[:60], result_canon[:60]))
    elif n_result < n_stored:
        fewer_slots.append((sym, n_stored, n_result, stored_canon[:60], result_canon[:60]))
    elif n_stored == n_result:
        same_slots_diff_numbering.append((sym, n_stored, stored_canon[:60], result_canon[:60]))

print(f"Extra slots (pipeline finds more):  {len(extra_slots)}")
print(f"Fewer slots (pipeline finds less):  {len(fewer_slots)}")
print(f"Same slot count, different struct:   {len(same_slots_diff_numbering)}")

if extra_slots:
    print(f"\n--- Extra slots ({len(extra_slots)}) ---")
    by_delta = {}
    for sym, ns, nr, s, r in extra_slots:
        d = nr - ns
        by_delta.setdefault(d, []).append(sym)
    for d in sorted(by_delta):
        syms = by_delta[d]
        print(f"  +{d} slots: {len(syms)} monomers")
        for s in syms[:5]:
            print(f"    {s}")
        if len(syms) > 5:
            print(f"    ... and {len(syms)-5} more")

if fewer_slots:
    print(f"\n--- Fewer slots ({len(fewer_slots)}) ---")
    for sym, ns, nr, s, r in fewer_slots:
        print(f"  {sym}: stored={ns} result={nr}")

if same_slots_diff_numbering:
    print(f"\n--- Same slot count ({len(same_slots_diff_numbering)}) ---")
    for sym, n, s, r in same_slots_diff_numbering[:20]:
        print(f"  {sym} ({n} slots)")
        print(f"    stored: {s}")
        print(f"    result: {r}")
