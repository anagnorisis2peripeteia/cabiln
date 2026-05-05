#!/usr/bin/env python3
"""Try to parse Leu_ol's mol block directly and find the issue."""
import pathlib, re
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = SDF.read_bytes().decode('utf-8')
records = text.split('$$$$')

# Find Leu_ol record index
leu_idx = None
for i, rec in enumerate(records):
    m = re.search(r'> +<m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m and m.group(1).strip() == 'Leu_ol':
        leu_idx = i
        print(f'Leu_ol is record #{i}')
        break

# Try to parse just that molblock
molblock = records[leu_idx].split('> <')[0].strip()
print(f'\nMolblock ({len(molblock)} chars):')
print(molblock)
print()

mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False, strictParsing=False)
print(f'MolFromMolBlock: {mol}')

# Try as a complete SDF record
single_sdf = records[leu_idx].strip() + '\n$$$$\n'
suppl = SDMolSupplier()
suppl.SetData(single_sdf.encode('utf-8'), removeHs=False, sanitize=False)
mols = [m for m in suppl]
print(f'SDMolSupplier from single record: {mols}')

# Check raw bytes around position of Leu_ol record
raw = SDF.read_bytes()
start_byte = text.encode('utf-8').index(records[leu_idx][:50].encode('utf-8'))
print(f'\nRaw bytes at record start (first 30): {raw[start_byte:start_byte+30]}')
