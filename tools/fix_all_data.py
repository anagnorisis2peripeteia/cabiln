#!/usr/bin/env python3
"""
Comprehensive SDF data fix script. Applies all corrections in order:
1. Data errors (gGlu, mePhe, C18FA, C16dAcid)
2. Chain length corrections
3. Fix dummy-on-O LG errors ([OH] → [H] where dummy bonds to O atom)
4. Safe slot additions (where pipeline finds extra slots, no structural change)
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

def parse_rgroups(rg):
    if isinstance(rg, list):
        return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]

def set_rgroups(sym, slot, new_lg):
    """Update a single LG value in the m_Rgroups string."""
    rg = df.loc[sym, 'm_Rgroups']
    parts = rg.split(',')
    while len(parts) < slot:
        parts.append('None')
    parts[slot - 1] = new_lg
    df.at[sym, 'm_Rgroups'] = ','.join(parts)

def canonical(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else None

LG_ATOMS = {
    '[H]': (1, 0), '[OH]': (8, 1),
    '[Cl]': (17, 0), '[Br]': (35, 0), '[I]': (53, 0),
}

SKIP_CAPS = {'Dmb', 'Ph', '_Bn', '_Et', '_Me', 'otbu', 'trt'}

# ═══════════════════════════════════════════════════════════════
# PHASE 1: Data error corrections
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PHASE 1: Data error corrections")
print("=" * 70)

CORRECTIONS = {
    'C18FA': '[2*]C(=O)CCCCCCCCCCCCCCCCCC(=O)O',
    'C16dAcid': '[2*]C(=O)CCCCCCCCCCCCCCCC(=O)O',
    'gGlu': '[1*]N([3*])[C@@H](C([4*])=O)CCC([2*])=O',
    'D_gGlu': '[1*]N([3*])[C@H](C([4*])=O)CCC([2*])=O',
    'mePhe_34diF': '[1*]N(C)[C@@H](Cc1ccc(F)c(F)c1)C([2*])=O',
    'D_mePhe_34diF': '[1*]N(C)[C@H](Cc1ccc(F)c(F)c1)C([2*])=O',
    'mePhe_3Cl': '[1*]N(C)[C@@H](Cc1cccc(Cl)c1)C([2*])=O',
    'D_mePhe_3Cl': '[1*]N(C)[C@H](Cc1cccc(Cl)c1)C([2*])=O',
}

for sym, new_smi in CORRECTIONS.items():
    new_mol = Chem.MolFromSmiles(new_smi)
    AllChem.Compute2DCoords(new_mol)
    df.at[sym, 'ROMol'] = new_mol
    print(f"  FIXED {sym}")

# ═══════════════════════════════════════════════════════════════
# PHASE 2: Chain length corrections
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("PHASE 2: Chain length corrections")
print("=" * 70)

EXPECTED = {
    'C14FA': 14, 'C16FA': 16, 'C18FA': 18, 'C20FA': 20, 'C22FA': 22,
    'C14dAcid': 14, 'C16dAcid': 16, 'C22dAcid': 22,
    'Dod_dAcid': 12,
    'Myr': 14, 'Pal': 16, 'Ste': 18,
}

chain_fixes = 0
for sym, expected_c in sorted(EXPECTED.items()):
    if sym not in df.index:
        continue
    mol = df.loc[sym, 'ROMol']
    actual_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    if actual_c == expected_c:
        continue

    smi = Chem.MolToSmiles(mol)
    is_diacid = sym not in ('Myr', 'Pal', 'Ste')

    if is_diacid:
        n_ch2 = expected_c - 2
        new_smi = '[2*]C(=O)' + 'C' * n_ch2 + 'C(=O)O'
    else:
        n_ch2 = expected_c - 2
        new_smi = '[2*]C(=O)' + 'C' * n_ch2 + 'C'

    new_mol = Chem.MolFromSmiles(new_smi)
    new_c = sum(1 for a in new_mol.GetAtoms() if a.GetAtomicNum() == 6)
    assert new_c == expected_c, f"{sym}: got {new_c}C, expected {expected_c}C"

    AllChem.Compute2DCoords(new_mol)
    df.at[sym, 'ROMol'] = new_mol
    chain_fixes += 1
    print(f"  FIXED {sym}: {actual_c}C -> {expected_c}C")

print(f"  {chain_fixes} chain lengths fixed")

# ═══════════════════════════════════════════════════════════════
# PHASE 3: Fix dummy-on-O LG errors
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("PHASE 3: Fix dummy-on-O LG errors")
print("=" * 70)

lg_fixes = 0
for sym in sorted(df.index):
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue
    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso < 1 or iso > len(rgroups):
            continue
        lg = rgroups[iso - 1]
        if lg != '[OH]':
            continue
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                set_rgroups(sym, iso, '[H]')
                lg_fixes += 1
                print(f"  FIXED {sym} slot {iso}: LG [OH] -> [H] (dummy bonds to O)")
                break

print(f"  {lg_fixes} LG errors fixed")

# ═══════════════════════════════════════════════════════════════
# PHASE 4: Safe slot additions
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("PHASE 4: Safe slot additions")
print("=" * 70)

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

def get_slot_numbers(smiles):
    return sorted(set(int(m) for m in re.findall(r'\[(\d+)\*\]', smiles)))

def strip_extra_slots(result_chuckles, stored_slots):
    mol = Chem.MolFromSmiles(result_chuckles)
    if mol is None:
        return None
    emol = Chem.RWMol(mol)
    to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            iso = atom.GetIsotope()
            if iso not in stored_slots:
                for nb in atom.GetNeighbors():
                    nb.SetNoImplicit(False)
                to_remove.append(atom.GetIdx())
    for idx in sorted(to_remove, reverse=True):
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except Exception:
        return None

safe_updates = 0
not_safe = 0

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

    stored_slots = set(get_slot_numbers(stored_canon))
    result_slots = set(get_slot_numbers(result_canon))
    extra_slots = result_slots - stored_slots

    if not extra_slots or not stored_slots.issubset(result_slots):
        not_safe += 1
        continue

    lg_ok = True
    for slot, lg in stored_lgs.items():
        if slot in result.leaving and result.leaving[slot] != lg:
            lg_ok = False
            break
    if not lg_ok:
        not_safe += 1
        continue

    stripped = strip_extra_slots(result.chuckles, stored_slots)
    if stripped != stored_canon:
        not_safe += 1
        continue

    new_mol = Chem.MolFromSmiles(result.chuckles)
    if new_mol is None:
        continue
    AllChem.Compute2DCoords(new_mol)

    max_slot = max(result.leaving.keys()) if result.leaving else 0
    max_slot = max(max_slot, 6)
    rg_parts = [result.leaving.get(i, 'None') for i in range(1, max_slot + 1)]
    df.at[sym, 'ROMol'] = new_mol
    df.at[sym, 'm_Rgroups'] = ','.join(rg_parts)
    safe_updates += 1

print(f"  {safe_updates} monomers safely updated with new slots")
print(f"  {not_safe} monomers skipped (structural differences)")

# ═══════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("SAVING")
print("=" * 70)
df_save = df.copy().reset_index()
df_save = df_save.rename(columns={'index': 'symbol'})
PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print("  SDF saved")
