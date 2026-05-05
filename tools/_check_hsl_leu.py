#!/usr/bin/env python3
"""Check status of Hsl and Leu_ol specifically."""
import pathlib
from rdkit import Chem
from rdkit.Chem import PandasTools

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
df = PandasTools.LoadSDF(str(SDF), removeHs=False, sanitize=False)
print(f'Total loaded: {len(df)}')
if 'm_abbr' in df.columns:
    for target in ['Hsl', 'Leu_ol', 'Lys_ol', 'Gly_ol']:
        rows = df[df['m_abbr'] == target]
        print(f'{target}: {len(rows)} rows')
        if len(rows) > 0:
            ct = rows.iloc[0].get('m_chem_types', 'MISSING')
            print(f'  m_chem_types: {ct}')

# Also check with ForwardSDMolSupplier
print('\nForwardSDMolSupplier:')
suppl = Chem.ForwardSDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=False)
results = {}
for mol in suppl:
    if mol is not None and mol.HasProp('m_abbr'):
        abbr = mol.GetProp('m_abbr')
        results[abbr] = results.get(abbr, 0) + 1

for target in ['Hsl', 'Leu_ol', 'Lys_ol']:
    print(f'{target}: {results.get(target, 0)} times')
print(f'Total unique: {len(results)}, total: {sum(results.values())}')
