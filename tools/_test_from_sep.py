#!/usr/bin/env python3
"""Test Leu_ol record starting from byte AFTER $$$$\\r\\n (simulating file-level parse)."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()

sep = b'$$$$'
leu_pos = raw.find(b'Leu_ol')
preceding = raw.rfind(sep, max(0, leu_pos - 2000), leu_pos)
next_sep = raw.find(sep, leu_pos)

print(f'Leu_ol at {leu_pos}, preceding sep at {preceding}, next sep at {next_sep}')

# What SDMolSupplier sees after consuming $$$$\r\n (6 bytes: 4 for $$$$, 2 for \r\n)
after_sep = preceding + len(sep) + 2  # skip $$$$\r\n
record_from_sep = raw[after_sep:next_sep + len(sep) + 1]  # include next $$$$\n
print(f'Record from after-sep: {record_from_sep[:60]!r}')

# Test with SetData
suppl = Chem.SDMolSupplier()
suppl.SetData(record_from_sep, removeHs=False, sanitize=False)
mols = list(suppl)
print(f'Parse from after-sep bytes: {mols}')
if mols and mols[0]:
    print(f'm_abbr: {mols[0].GetProp("m_abbr")}')

# Also test Hsl and Lys_ol the same way for comparison
text = raw.decode('utf-8')
records = text.split('$$$$')
for target in ['Hsl', 'Leu_ol', 'Lys_ol']:
    import re
    for i, rec in enumerate(records):
        m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
        if m and m.group(1).strip() == target:
            # Calculate byte position of preceding $$$$
            # Each record in the text corresponds to one SDF record
            # We need to find the byte position in raw
            preceding_text = '$$$$'.join(records[:i])
            sep_byte_pos = len(preceding_text.encode('utf-8')) + (i * 4)  # estimate
            # Better: search in raw
            search_pos = raw.find(target.encode('utf-8'))
            prec_sep = raw.rfind(b'$$$$', max(0, search_pos - 2000), search_pos)
            next_s = raw.find(b'$$$$', search_pos)
            after = prec_sep + 6  # skip $$$$\r\n
            rec_bytes = raw[after:next_s + 5]
            s2 = Chem.SDMolSupplier()
            s2.SetData(rec_bytes, removeHs=False, sanitize=False)
            ml = list(s2)
            ok = ml and ml[0] is not None
            print(f'{target} (record {i}): parse_ok={ok}, start={rec_bytes[:30]!r}')
            break
