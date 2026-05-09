#!/usr/bin/env python3
"""
Add TATA and TBAB scaffold crosslinkers to monomers.csv + monomers.sdf.

  TATA — 1,3,5-Triacryloyl-1,3,5-triazinane
    3 Cys react via thio-Michael addition at the beta-carbons of acrylamide arms.
    No leaving group (addition chemistry).

  TBAB — N,N',N''-Benzene-1,3,5-triyl-tris(2-bromoacetamide)
    3 Cys react via thioether alkylation at CH2-Br arms (same as TBMB mechanism,
    but benzene NH-amide arms instead of direct CH2).
    Leaving group: [Br].
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

CSV = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.csv')
SDF = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.sdf')

NEW_MONOMERS = [
    {
        'token':     'TATA',
        'input':     'O=C(C=C)N1CN(C(=O)C=C)CN(C(=O)C=C)C1',
        'name':      '1,3,5-Triacryloyl-1,3,5-triazinane',
        'type':      'linker',
        'synonyms':  '',
        'chuckles':  '[4*]CCC(=O)N1CN(C(=O)CC[5*])CN(C(=O)CC[6*])C1',
        'r1_leaving': '', 'r2_leaving': '', 'r3_leaving': '',
        'r4_leaving': '', 'r5_leaving': '', 'r6_leaving': '',
        'm_subtype':   'crosslink',
        'm_chem_types': '4:thia_michael_c,5:thia_michael_c,6:thia_michael_c',
    },
    {
        'token':     'TBAB',
        'input':     'BrCC(=O)Nc1cc(NC(=O)CBr)cc(NC(=O)CBr)c1',
        'name':      "N,N',N''-Benzene-1,3,5-triyl-tris(2-bromoacetamide)",
        'type':      'linker',
        'synonyms':  '',
        'chuckles':  '[4*]CC(=O)Nc1cc(NC(=O)C[5*])cc(NC(=O)C[6*])c1',
        'r1_leaving': '', 'r2_leaving': '', 'r3_leaving': '',
        'r4_leaving': '[Br]', 'r5_leaving': '[Br]', 'r6_leaving': '[Br]',
        'm_subtype':   'crosslink',
        'm_chem_types': '4:alkyl_halide_c,5:alkyl_halide_c,6:alkyl_halide_c',
    },
]

# ---------------------------------------------------------------------------
# 1. Patch CSV
# ---------------------------------------------------------------------------
print('Patching CSV...')
df = pd.read_csv(CSV)

for mon in NEW_MONOMERS:
    tok = mon['token']
    if tok in df['token'].values:
        print(f'  {tok}: already in CSV — skipping')
        continue
    row = {
        'token':      tok,
        'input':      mon['input'],
        'name':       mon['name'],
        'type':       mon['type'],
        'synonyms':   mon['synonyms'],
        'chuckles':   mon['chuckles'],
        'r1_leaving': mon['r1_leaving'],
        'r2_leaving': mon['r2_leaving'],
        'r3_leaving': mon['r3_leaving'],
        'r4_leaving': mon['r4_leaving'],
        'r5_leaving': mon['r5_leaving'],
        'r6_leaving': mon['r6_leaving'],
    }
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    print(f'  {tok}: added to CSV')

df.to_csv(CSV, index=False)
print('CSV saved.\n')

# ---------------------------------------------------------------------------
# 2. Generate SDF mol blocks from CHUCKLES and append
# ---------------------------------------------------------------------------
print('Patching SDF...')

with open(SDF, 'r', encoding='utf-8') as f:
    sdf_content = f.read()

# Check which tokens already exist in SDF
existing_abbrs = set(re.findall(r'>  <m_abbr>[^\n]*\n(.+?)\n', sdf_content))
print(f'  SDF has {sdf_content.count("$$$$")} existing records')

new_records = []
for mon in NEW_MONOMERS:
    tok = mon['token']
    if tok in existing_abbrs:
        print(f'  {tok}: already in SDF — skipping')
        continue

    mol = Chem.MolFromSmiles(mon['chuckles'])
    if mol is None:
        print(f'  {tok}: ERROR — invalid CHUCKLES: {mon["chuckles"]}')
        continue

    AllChem.Compute2DCoords(mol)

    # Build m_Rgroups string: 6 slots, None for empty
    rg_list = [
        mon['r1_leaving'] or 'None',
        mon['r2_leaving'] or 'None',
        mon['r3_leaving'] or 'None',
        mon['r4_leaving'] or 'None',
        mon['r5_leaving'] or 'None',
        mon['r6_leaving'] or 'None',
    ]
    m_rgroups = ', '.join(rg_list)

    # Get next record number from SDF
    rec_num = sdf_content.count('$$$$') + len(new_records) + 1

    mol_block = Chem.MolToMolBlock(mol)
    # Set name line to token
    lines = mol_block.split('\n')
    lines[0] = tok
    mol_block = '\n'.join(lines)

    props = (
        f'>  <symbol>  ({rec_num}) \n{tok}\n\n'
        f'>  <m_Rgroups>  ({rec_num}) \n{m_rgroups}\n\n'
        f'>  <m_abbr>  ({rec_num}) \n{tok}\n\n'
        f'>  <m_name>  ({rec_num}) \n{mon["name"]}\n\n'
        f'>  <m_type>  ({rec_num}) \n{mon["type"]}\n\n'
        f'>  <m_subtype>  ({rec_num}) \n{mon["m_subtype"]}\n\n'
        f'>  <m_chem_types>  ({rec_num}) \n{mon["m_chem_types"]}\n\n'
    )
    record = mol_block + props + '$$$$\n'
    new_records.append((tok, record))
    print(f'  {tok}: built mol block ({mol.GetNumAtoms()} atoms)')

if new_records:
    with open(SDF, 'a', encoding='utf-8', newline='\n') as f:
        for tok, rec in new_records:
            f.write(rec)
    print(f'\nSDF updated — appended {len(new_records)} records')
else:
    print('\nSDF unchanged — all monomers already present')

# ---------------------------------------------------------------------------
# 3. Verify
# ---------------------------------------------------------------------------
print('\n=== Verification ===\n')
df2 = pd.read_csv(CSV)
suppl = Chem.SDMolSupplier(SDF, removeHs=False)
sdf_abbrs = {}
for m in suppl:
    if m is None: continue
    abbr = m.GetPropsAsDict().get('m_abbr', '').strip()
    if abbr:
        sdf_abbrs[abbr] = m

for mon in NEW_MONOMERS:
    tok = mon['token']
    csv_row = df2[df2['token'] == tok]
    in_csv = not csv_row.empty
    in_sdf = tok in sdf_abbrs

    print(f'{tok}:')
    print(f'  CSV: {in_csv}  chuckles={csv_row["chuckles"].iloc[0] if in_csv else "MISSING"}')
    print(f'  SDF: {in_sdf}')
    if in_sdf:
        sdf_mol = sdf_abbrs[tok]
        dummies = sorted(a.GetIsotope() for a in sdf_mol.GetAtoms() if a.GetAtomicNum() == 0)
        chem_types = sdf_mol.GetPropsAsDict().get('m_chem_types', '')
        print(f'  dummy R-groups in mol: {dummies}')
        print(f'  m_chem_types: {chem_types}')

    # Scaffold detection check
    from tools.live_renderer import _s2c_get_scaffold_patterns, _sdf_cache
    _sdf_cache['mols'] = None  # force reload
    patterns = _s2c_get_scaffold_patterns.__wrapped__() if hasattr(_s2c_get_scaffold_patterns, '__wrapped__') else None
    print()
