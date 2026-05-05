#!/usr/bin/env python3
"""Restore Leu_ol record start to match normal record pattern (two CRLFs after $$$$).
Also diagnose the actual failure by comparing with neighboring records."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

# Find Leu_ol
leu_idx = None
for i, rec in enumerate(records):
    if 'Leu_ol' in rec and '>  <m_abbr>' in rec:
        leu_idx = i
        break

leu_rec = records[leu_idx]
print(f'Leu_ol index: {leu_idx}')
print(f'Current start: {leu_rec[:40]!r}')

# Restore to standard pattern: \r\n\r\n     RDKit
if leu_rec.startswith('\r\n     RDKit'):
    # Currently single CRLF (my broken fix) - restore to double CRLF
    records[leu_idx] = '\r\n' + leu_rec
    print('Restored to double CRLF pattern')
elif leu_rec.startswith('\r\n\r\n     RDKit'):
    print('Already double CRLF - no change needed')
else:
    print(f'Unexpected start: {leu_rec[:20]!r}')

# Write back
SDF.write_bytes('$$$$'.join(records).encode('utf-8'))
print('Written back.')

# Now test: what does a working neighbor record look like?
print('\n--- Comparing records ---')
for i in [leu_idx - 1, leu_idx, leu_idx + 1]:
    rec = records[i]
    # Find m_abbr
    import re
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    abbr = m.group(1).strip() if m else '?'
    start = rec[:50]

    # Test isolated parse
    test_sdf = (rec + '$$$$\n').encode('utf-8')
    suppl = Chem.SDMolSupplier()
    suppl.SetData(test_sdf, removeHs=False, sanitize=False)
    mols = list(suppl)
    ok = mols and mols[0] is not None
    print(f'Record {i} ({abbr}): start={start!r}, parse_ok={ok}')
