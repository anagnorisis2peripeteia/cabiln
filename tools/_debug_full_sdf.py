#!/usr/bin/env python3
"""Debug full SDF loading - find exactly which record fails and why."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

# Use ForwardSDMolSupplier which processes one at a time
suppl = Chem.ForwardSDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=False)

count = 0
none_count = 0
leu_found = False
for i, mol in enumerate(suppl):
    count += 1
    if mol is None:
        none_count += 1
        print(f'Record {i}: None')
    else:
        abbr = mol.GetProp('m_abbr') if mol.HasProp('m_abbr') else '?'
        if abbr == 'Leu_ol':
            leu_found = True
            print(f'Record {i}: {abbr} <- FOUND!')
        if i >= 308 and i <= 316:
            print(f'Record {i}: {abbr}')

print(f'\nTotal: {count} records, {none_count} None, leu_found={leu_found}')
