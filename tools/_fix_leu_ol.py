#!/usr/bin/env python3
"""Fix Leu_ol's double-blank-line header by removing the extra CRLF at record start.

The issue: after splitting on $$$$, Leu_ol's record starts with \\r\\n\\r\\n     RDKit
(two CRLFs before the mol block header). V2000 needs exactly 3 header lines before the
counts line, so two leading blank lines push the counts to line 5, causing parse failure.

Fix: remove one \\r\\n from the start of Leu_ol's record in the SDF, so the record
starts with \\r\\n     RDKit (one blank = mol name, then RDKit = program, blank = comment,
counts = line 4).
"""
import pathlib, re

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

# Find Leu_ol record and check its start
leu_idx = None
for i, rec in enumerate(records):
    if 'Leu_ol' in rec and '>  <m_abbr>' in rec:
        leu_idx = i
        break

if leu_idx is None:
    print('ERROR: Leu_ol not found')
    raise SystemExit(1)

leu_rec = records[leu_idx]
print(f'Leu_ol record index: {leu_idx}')
print(f'Record starts with: {leu_rec[:30]!r}')

# Check for double CRLF at start
if leu_rec.startswith('\r\n\r\n'):
    print('Found double CRLF at start - fixing...')
    records[leu_idx] = leu_rec[2:]  # Remove one \r\n (2 bytes)
    fixed_text = '$$$$'.join(records)
    SDF.write_bytes(fixed_text.encode('utf-8'))
    print('Fixed! Wrote updated SDF.')
elif leu_rec.startswith('\r\n'):
    print('Single CRLF at start - already correct, no fix needed.')
else:
    print(f'Unexpected start: {leu_rec[:20]!r}')

# Verify the fix
raw2 = SDF.read_bytes()
text2 = raw2.decode('utf-8')
records2 = text2.split('$$$$')
leu_rec2 = records2[leu_idx]
print(f'After fix, record starts with: {leu_rec2[:40]!r}')

# Quick parse check
from rdkit import Chem
suppl = Chem.SDMolSupplier()
test_sdf = leu_rec2 + '$$$$\n'
suppl.SetData(test_sdf.encode('utf-8'), removeHs=False, sanitize=False)
mols = list(suppl)
print(f'Parse check: {mols}')
if mols and mols[0] is not None:
    props = mols[0].GetPropsAsDict()
    print(f'm_abbr: {props.get("m_abbr", "NOT FOUND")}')
