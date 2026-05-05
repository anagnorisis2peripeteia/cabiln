#!/usr/bin/env python3
"""Find all records that fail to parse and show their first few lines."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

# Use SDMolSupplier with strict parsing to identify problem records
suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=True)
failed_indices = []
for i, mol in enumerate(suppl):
    if mol is None:
        failed_indices.append(i)

print(f'Failed record indices (0-based): {failed_indices}')
print(f'Total records attempted: {i+1}')

# Show context around failed records using text
text = SDF.read_bytes().decode('utf-8')
records = text.split('$$$$')

import re
def gp(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else '?'

for idx in failed_indices:
    rec = records[idx] if idx < len(records) else '(out of range)'
    abbr = gp(rec, 'm_abbr')
    # Show first 6 lines of mol block
    lines = rec.strip().split('\n')[:8]
    print(f'\nRecord {idx} ({abbr}):')
    for l in lines:
        print(f'  {repr(l)}')
