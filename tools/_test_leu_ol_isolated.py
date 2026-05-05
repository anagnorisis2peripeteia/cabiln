#!/usr/bin/env python3
"""Extract Leu_ol as standalone SDF and test loading in various ways."""
import pathlib, tempfile, os
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

leu_rec = None
leu_idx = None
for i, rec in enumerate(records):
    if 'Leu_ol' in rec and '>  <m_abbr>' in rec:
        leu_rec = rec
        leu_idx = i
        print(f'Leu_ol is record index {i}')
        break

# Write Leu_ol as a standalone SDF, exactly as it would be read from file
# (with leading \r\n intact)
leu_sdf = leu_rec + '$$$$\n'
leu_sdf_bytes = leu_sdf.encode('utf-8')
print(f'SDF bytes: {leu_sdf_bytes[:100]!r}')

with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False, mode='wb') as f:
    f.write(leu_sdf_bytes)
    tmppath = f.name

print(f'\nTesting SDMolSupplier on isolated SDF: {tmppath}')
suppl = Chem.SDMolSupplier(tmppath, removeHs=False, sanitize=False, strictParsing=False)
mols = list(suppl)
print(f'  strictParsing=False: {mols}')

# Also try with explicit CRLF bytes
leu_sdf_crlf = leu_rec.encode('utf-8') + b'$$$$\r\n'
with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False, mode='wb') as f:
    f.write(leu_sdf_crlf)
    tmppath2 = f.name

print(f'\nTesting with CRLF terminator: {tmppath2}')
suppl2 = Chem.SDMolSupplier(tmppath2, removeHs=False, sanitize=False, strictParsing=False)
mols2 = list(suppl2)
print(f'  Result: {mols2}')

# Try with SetData
print('\nTesting SDMolSupplier.SetData:')
suppl3 = Chem.SDMolSupplier()
suppl3.SetData(leu_sdf_crlf, removeHs=False, sanitize=False)
mols3 = list(suppl3)
print(f'  Result: {mols3}')
if mols3 and mols3[0] is not None:
    props = list(mols3[0].GetPropsAsDict().keys())
    print(f'  Properties loaded: {props}')

# Try a fixed version: strip leading \r\n and add back 3-line header
print('\nTesting with corrected header (3 lines before counts):')
rec_stripped = leu_rec.lstrip('\r\n ')
# Ensure 3-line header: name\nprogram\ncomment\ncounts
fixed_rec = '\n' + rec_stripped  # single blank line = mol name, then RDKit header = program, then blank = comment
fixed_sdf = fixed_rec + '$$$$\n'
suppl4 = Chem.SDMolSupplier()
suppl4.SetData(fixed_sdf.encode('utf-8'), removeHs=False, sanitize=False)
mols4 = list(suppl4)
print(f'  Result: {mols4}')

os.unlink(tmppath)
os.unlink(tmppath2)
