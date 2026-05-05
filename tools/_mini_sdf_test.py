#!/usr/bin/env python3
"""Extract records 309-314 as mini SDF and test loading to isolate the failure."""
import pathlib, tempfile, os
from rdkit import Chem
from rdkit.Chem import PandasTools

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

# Find Leu_ol index
import re
leu_idx = next(i for i, r in enumerate(records)
               if re.search(r'>  <m_abbr>.*?\r?\n *Leu_ol', r))
print(f'Leu_ol at index {leu_idx}')

# Show neighbors
for i in range(leu_idx - 3, leu_idx + 3):
    if 0 <= i < len(records):
        m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', records[i])
        abbr = m.group(1).strip() if m else '?'
        print(f'  Record {i}: {abbr}')

# Build mini SDF with records 309-314
start = leu_idx - 3
end = leu_idx + 3
mini_records = records[start:end+1]

# The first record needs special handling - it shouldn't have a leading \r\n\r\n
# because there's no preceding $$$$ in the mini SDF
# BUT the records all start with \r\n\r\n which confuses SetData...
# Use the proper file representation: first record gets its leading \r\n from after $$$$

# Write as proper SDF file
mini_sdf_text = '$$$$'.join(mini_records) + '$$$$\n'
mini_sdf_bytes = mini_sdf_text.encode('utf-8')

with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False, mode='wb') as f:
    f.write(mini_sdf_bytes)
    tmppath = f.name

print(f'\nMini SDF: {tmppath}')
print(f'Mini SDF starts with: {mini_sdf_bytes[:60]!r}')

# Test with SDMolSupplier
suppl = Chem.SDMolSupplier(tmppath, removeHs=False, sanitize=False, strictParsing=False)
mols = list(suppl)
print(f'Records loaded: {len(mols)}')
abbrs = []
for mol in mols:
    if mol is not None:
        abbrs.append(mol.GetProp('m_abbr') if mol.HasProp('m_abbr') else '?')
    else:
        abbrs.append('None')
print(f'Abbrs: {abbrs}')

# Test with PandasTools
df = PandasTools.LoadSDF(tmppath, removeHs=False, sanitize=False)
print(f'PandasTools loaded: {len(df)} records')
if 'm_abbr' in df.columns:
    print(f'Abbrs: {df["m_abbr"].tolist()}')
    print(f'Leu_ol present: {"Leu_ol" in df["m_abbr"].values}')

os.unlink(tmppath)
