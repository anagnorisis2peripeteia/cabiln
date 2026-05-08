#!/usr/bin/env python3
"""
Fix three genuine stereo errors:
  meT   – NMe-L-Thr: CHUCKLES assembles to D-form (alpha=R); must be L (alpha=S)
  D_meT – NMe-D-Thr: CHUCKLES assembles to L-form (alpha=S); must be D (alpha=R)
  D_Aze – D-Aze-2-CA: identical CHUCKLES/input to L-Aze; must be D-form (alpha=R)

Updates both monomers.csv and monomers.sdf.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem

CSV = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.csv')
SDF = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.sdf')

# ---------------------------------------------------------------------------
# Verified correct values (empirically derived in _fix_stereo_errors.py tests)
# ---------------------------------------------------------------------------
FIXES = {
    'meT': {
        'input':    'CN[C@@H]([C@H](O)C)C(=O)O',   # NMe-L-Thr (2S,3R)
        'chuckles': '[1*]N(C)[C@H](C([2*])=O)[C@@H](C)O[3*]',
    },
    'D_meT': {
        'input':    'CN[C@H]([C@@H](O)C)C(=O)O',   # NMe-D-Thr (2R,3S)
        'chuckles': '[1*]N(C)[C@@H](C([2*])=O)[C@H](C)O[3*]',
    },
    'D_Aze': {
        'input':    'O=C(O)[C@H]1CCN1',              # D-Aze-2-CA (R)
        'chuckles': '[1*]N1CC[C@@H]1C([2*])=O',
    },
}

LG_ATOMS = {'[H]': (1, 0), '[OH]': (8, 1), '[Cl]': (17, 0), '[Br]': (35, 0), '[I]': (53, 0)}

def restore_full_smiles(mol, rgroups):
    emol = Chem.RWMol(mol)
    try:
        Chem.Kekulize(emol, clearAromaticFlags=True)
    except Exception:
        pass
    dummies = [(a.GetIdx(), a.GetIsotope()) for a in emol.GetAtoms() if a.GetAtomicNum() == 0]
    remove = []
    for dummy_idx, slot in dummies:
        if slot < 1 or slot > len(rgroups):
            continue
        lg = rgroups[slot - 1]
        if lg is None:
            continue
        if lg == '[H]':
            for nb in emol.GetAtomWithIdx(dummy_idx).GetNeighbors():
                nb.SetNoImplicit(False)
            remove.append(dummy_idx)
        elif lg == '[OH]':
            a = emol.GetAtomWithIdx(dummy_idx)
            a.SetAtomicNum(8); a.SetIsotope(0); a.SetNumExplicitHs(1)
        elif lg in ('[Cl]', '[Br]', '[I]'):
            a = emol.GetAtomWithIdx(dummy_idx)
            a.SetAtomicNum(LG_ATOMS[lg][0]); a.SetIsotope(0)
        else:
            remove.append(dummy_idx)
    for idx in sorted(remove, reverse=True):
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except Exception:
        return None

def cip_dict(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return {}
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return {a.GetIdx(): a.GetPropsAsDict()['_CIPCode'] for a in mol.GetAtoms() if a.HasProp('_CIPCode')}

# ---------------------------------------------------------------------------
# Verify before applying
# ---------------------------------------------------------------------------
print('=== Pre-fix verification ===')
df = pd.read_csv(CSV).set_index('token')
rgroups_map = {}
for tok in FIXES:
    row = df.loc[tok]
    rgroups_map[tok] = [row.get(f'r{i}_leaving') for i in range(1, 7)
                        if isinstance(row.get(f'r{i}_leaving'), str)]

print('\nAssembled SMILES from NEW CHUCKLES:')
for tok, fix in FIXES.items():
    mol = Chem.MolFromSmiles(fix['chuckles'])
    rg = rgroups_map.get(tok, ['[H]', '[OH]'])
    asm = restore_full_smiles(mol, rg) if mol else None
    cips = cip_dict(asm) if asm else {}
    print(f'  {tok:8s}  {asm}  CIPs={cips}')
    input_cips = cip_dict(fix['input'])
    print(f'           input {fix["input"]}  CIPs={input_cips}')
print()

# Abort if CHUCKLES and input don't produce the same connectivity (flat)
errors = []
for tok, fix in FIXES.items():
    mol = Chem.MolFromSmiles(fix['chuckles'])
    rg = rgroups_map.get(tok, ['[H]', '[OH]'])
    asm = restore_full_smiles(mol, rg) if mol else None
    if asm is None:
        errors.append(f'{tok}: restore failed')
        continue
    # Check flat structure matches input
    def flat(s):
        m = Chem.MolFromSmiles(s)
        if not m: return None
        Chem.RemoveStereochemistry(m)
        return Chem.MolToSmiles(m)
    if flat(asm) != flat(fix['input']):
        errors.append(f'{tok}: assembled flat != input flat\n  asm={asm}\n  input={fix["input"]}')

if errors:
    print('ERRORS – aborting:')
    for e in errors:
        print(' ', e)
    sys.exit(1)

print('All pre-checks passed. Applying fixes...\n')

# ---------------------------------------------------------------------------
# Update CSV
# ---------------------------------------------------------------------------
rows = []
with open(CSV, newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows.append(header)
    col_input    = header.index('input')
    col_chuckles = header.index('chuckles')
    for row in reader:
        tok = row[0]
        if tok in FIXES:
            old_input    = row[col_input]
            old_chuckles = row[col_chuckles]
            row[col_input]    = FIXES[tok]['input']
            row[col_chuckles] = FIXES[tok]['chuckles']
            print(f'CSV {tok}:')
            print(f'  input    {old_input} -> {row[col_input]}')
            print(f'  chuckles {old_chuckles} -> {row[col_chuckles]}')
        rows.append(row)

with open(CSV, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerows(rows)
print('\nCSV updated.')

# ---------------------------------------------------------------------------
# Update SDF
# ---------------------------------------------------------------------------
print('\nLoading SDF...')
sdf_df = PandasTools.LoadSDF(SDF)
sdf_df = sdf_df.set_index('m_abbr')
print(f'  {len(sdf_df)} entries')

for tok, fix in FIXES.items():
    if tok not in sdf_df.index:
        print(f'  {tok}: NOT IN SDF – skipping SDF update')
        continue
    new_mol = Chem.MolFromSmiles(fix['chuckles'])
    if new_mol is None:
        print(f'  {tok}: mol parse failed – skipping')
        continue
    AllChem.Compute2DCoords(new_mol)
    sdf_df.at[tok, 'ROMol'] = new_mol
    print(f'  SDF {tok}: updated mol block from {fix["chuckles"]}')

print('\nWriting SDF...')
PandasTools.WriteSDF(sdf_df.reset_index(), SDF,
                     molColName='ROMol', idName='m_abbr',
                     properties=list(sdf_df.columns.drop('ROMol')))
print('Done.')
