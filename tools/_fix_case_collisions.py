#!/usr/bin/env python3
"""
Fix three case-insensitive token collision pairs:
  1. Boc.input  — contains CHUCKLES notation instead of free SMILES
  2. Hser       — input/chuckles encode Serine (missing CH2), should be Homoserine
  3. otbu       — input/chuckles missing the ether oxygen
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

CSV = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.csv')
SDF = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.sdf')

# CSV fixes: token -> {field: new_value}
CSV_FIXES = {
    'Boc': {
        'input': 'CC(C)(C)OC(=O)O',   # free SMILES (same structure as boc.input)
    },
    'Hser': {
        'input':    'N[C@@H](CCO)C(=O)O',                   # Homoserine (was Serine)
        'chuckles': '[1*]N([3*])[C@@H](CCO[4*])C([2*])=O',  # Homoserine (was Serine)
    },
    'otbu': {
        'input':    'OC(C)(C)C',       # tert-butanol free form (was missing O)
        'chuckles': '[1*]OC(C)(C)C',   # O-tert-butyl ether (was missing O)
    },
    'OtBu': {
        'input': 'OC(C)(C)C',          # free SMILES (was CHUCKLES [1*]OC(C)(C)C)
    },
}

# SDF mol-block fixes: only tokens whose chuckles actually changed
SDF_CHUCKLES_FIXES = {
    'Hser': '[1*]N([3*])[C@@H](CCO[4*])C([2*])=O',
    'otbu': '[1*]OC(C)(C)C',
}

# ---------------------------------------------------------------------------
# 1. Patch CSV
# ---------------------------------------------------------------------------
print('Patching CSV...')
df = pd.read_csv(CSV)
for tok, fields in CSV_FIXES.items():
    mask = df['token'] == tok
    if not mask.any():
        print(f'  WARNING: {tok} not found in CSV')
        continue
    for field, val in fields.items():
        old = df.loc[mask, field].iloc[0]
        df.loc[mask, field] = val
        print(f'  {tok}.{field}: {old!r} -> {val!r}')

df.to_csv(CSV, index=False)
print('CSV saved.\n')

# ---------------------------------------------------------------------------
# 2. Patch SDF mol blocks (only structural changes)
# ---------------------------------------------------------------------------

def smiles_to_molblock(smi, name=''):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smi}')
    AllChem.Compute2DCoords(mol)
    mb = Chem.MolToMolBlock(mol, forceV3000=False)
    lines = mb.split('\n')
    lines[0] = name
    return '\n'.join(lines)

new_molblocks = {}
for tok, smi in SDF_CHUCKLES_FIXES.items():
    mb = smiles_to_molblock(smi, tok)
    new_molblocks[tok] = mb
    print(f'{tok}: built mol block from {smi}')

print(f'\nPatching SDF ({len(new_molblocks)} records to update)...')
with open(SDF, 'r', encoding='utf-8') as f:
    content = f.read()

records = content.split('$$$$\n')
if records[-1].strip() == '':
    records = records[:-1]

print(f'  Loaded {len(records)} SDF records')

updated = []
patched = set()
for rec in records:
    m = re.search(r'>  <m_abbr>[^\n]*\n(.+?)\n', rec)
    tok = m.group(1).strip() if m else None

    if tok in new_molblocks:
        end_idx = rec.find('M  END')
        if end_idx == -1:
            print(f'  WARNING: M  END not found for {tok}')
            updated.append(rec)
            continue
        prop_block = rec[end_idx + len('M  END'):]
        new_rec = new_molblocks[tok] + prop_block
        updated.append(new_rec)
        patched.add(tok)
        print(f'  Patched mol block for {tok}')
    else:
        updated.append(rec)

missing = set(new_molblocks) - patched
if missing:
    print(f'\nWARNING: tokens not found in SDF: {missing}')

with open(SDF, 'w', encoding='utf-8', newline='\n') as f:
    for rec in updated:
        if not rec.endswith('\n'):
            rec += '\n'
        f.write(rec + '$$$$\n')

print(f'SDF written with {len(updated)} records. Patched: {patched}\n')

# ---------------------------------------------------------------------------
# 3. Verify
# ---------------------------------------------------------------------------
print('=== Verification ===\n')
df2 = pd.read_csv(CSV)
for tok in CSV_FIXES:
    row = df2[df2['token'] == tok].iloc[0]
    print(f'{tok}:')
    print(f'  input    = {row["input"]}')
    print(f'  chuckles = {row["chuckles"]}')
    mol = Chem.MolFromSmiles(row['input'])
    print(f'  input valid: {mol is not None}')
    print()
