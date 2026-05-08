#!/usr/bin/env python3
"""
Surgically replace mol blocks in monomers.sdf for specific tokens.
Preserves all >  <property> fields unchanged.
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import AllChem

SDF = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.sdf')

# New CHUCKLES mol blocks (token -> SMILES)
NEW_CHUCKLES = {
    'meT':   '[1*]N(C)[C@H](C([2*])=O)[C@@H](C)O[3*]',
    'D_meT': '[1*]N(C)[C@@H](C([2*])=O)[C@H](C)O[3*]',
    'D_Aze': '[1*]N1CC[C@@H]1C([2*])=O',
}

def smiles_to_molblock(smi, name=''):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smi}')
    AllChem.Compute2DCoords(mol)
    mb = Chem.MolToMolBlock(mol, forceV3000=False)
    # Replace the default blank name line with the monomer abbreviation
    lines = mb.split('\n')
    lines[0] = name
    return '\n'.join(lines)

# Build lookup of new mol blocks
new_molblocks = {}
for tok, smi in NEW_CHUCKLES.items():
    mb = smiles_to_molblock(smi, tok)
    new_molblocks[tok] = mb
    print(f'{tok}: built mol block from {smi}')

# Parse SDF into records
with open(SDF, 'r', encoding='utf-8') as f:
    content = f.read()

records = content.split('$$$$\n')
# Last split might be empty
if records[-1].strip() == '':
    records = records[:-1]

print(f'\nLoaded {len(records)} SDF records')

updated = []
patched = set()
for rec in records:
    # Find the m_abbr property to identify the token
    m = re.search(r'>  <m_abbr>[^\n]*\n(.+?)\n', rec)
    tok = m.group(1).strip() if m else None

    if tok in new_molblocks:
        # Split at M  END to separate mol block from property block
        end_idx = rec.find('M  END')
        if end_idx == -1:
            print(f'WARNING: M  END not found in record for {tok}')
            updated.append(rec)
            continue
        prop_block = rec[end_idx + len('M  END'):]  # everything after M  END
        new_rec = new_molblocks[tok] + prop_block
        updated.append(new_rec)
        patched.add(tok)
        print(f'Patched mol block for {tok}')
    else:
        updated.append(rec)

missing = set(new_molblocks) - patched
if missing:
    print(f'\nWARNING: tokens not found in SDF: {missing}')

# Write back
with open(SDF, 'w', encoding='utf-8', newline='\n') as f:
    for rec in updated:
        if not rec.endswith('\n'):
            rec += '\n'
        f.write(rec + '$$$$\n')

print(f'\nSDF written with {len(updated)} records. Patched: {patched}')

# Verify by loading back
from rdkit.Chem import PandasTools
df = PandasTools.LoadSDF(SDF)
print(f'Reload check: {len(df)} entries, columns={df.columns.tolist()[:5]}')
for tok in new_molblocks:
    rows = df[df['m_abbr'] == tok]
    if rows.empty:
        print(f'  {tok}: NOT FOUND after reload')
    else:
        mol = rows.iloc[0]['ROMol']
        smi = Chem.MolToSmiles(mol) if mol else None
        print(f'  {tok}: reloaded SMILES = {smi}')
