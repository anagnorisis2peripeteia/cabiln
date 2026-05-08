#!/usr/bin/env python3
"""
Patch T (L-Thr) and DThr (D-Thr) SDF mol blocks.
The SDF has wrong beta-carbon stereo; CSV is correct.
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import AllChem

SDF = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.sdf')

# CSV-correct CHUCKLES (verified: T=2S,3R; DThr=2R,3S)
NEW_CHUCKLES = {
    'T':    '[1*]N([3*])[C@H](C([2*])=O)[C@@H](C)O[4*]',
    'DThr': '[1*]N([3*])[C@@H](C([2*])=O)[C@H](C)O[4*]',
}

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
for tok, smi in NEW_CHUCKLES.items():
    mb = smiles_to_molblock(smi, tok)
    new_molblocks[tok] = mb
    print(f'{tok}: built from {smi}')

with open(SDF, 'r', encoding='utf-8') as f:
    content = f.read()

records = content.split('$$$$\n')
if records[-1].strip() == '':
    records = records[:-1]

print(f'\nLoaded {len(records)} SDF records')

updated = []
patched = set()
for rec in records:
    m = re.search(r'>  <m_abbr>[^\n]*\n(.+?)\n', rec)
    tok = m.group(1).strip() if m else None
    if tok in new_molblocks:
        end_idx = rec.find('M  END')
        if end_idx == -1:
            print(f'WARNING: M  END not found for {tok}')
            updated.append(rec)
            continue
        prop_block = rec[end_idx + len('M  END'):]
        updated.append(new_molblocks[tok] + prop_block)
        patched.add(tok)
        print(f'Patched {tok}')
    else:
        updated.append(rec)

missing = set(new_molblocks) - patched
if missing:
    print(f'WARNING: tokens not found in SDF: {missing}')

with open(SDF, 'w', encoding='utf-8', newline='\n') as f:
    for rec in updated:
        if not rec.endswith('\n'):
            rec += '\n'
        f.write(rec + '$$$$\n')

print(f'\nSDF written. Patched: {patched}')

# Verify
suppl = Chem.SDMolSupplier(SDF, removeHs=False)
for mol in suppl:
    if mol is None:
        continue
    abbr = mol.GetPropsAsDict().get('m_abbr', '').strip() if mol.HasProp('m_abbr') else ''
    if abbr not in new_molblocks:
        continue
    sdf_smi = Chem.MolToSmiles(mol)
    exp_smi = Chem.MolToSmiles(Chem.MolFromSmiles(NEW_CHUCKLES[abbr]))
    match = sdf_smi == exp_smi
    print(f'{abbr}: SDF={sdf_smi}  match={match}')
