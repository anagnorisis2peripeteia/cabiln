#!/usr/bin/env python3
"""Check actual file bytes around Leu_ol's record start."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()

# Find the preceding $$$$ in the file bytes
leu_text = b'Leu_ol'
leu_pos = raw.find(leu_text)
print(f'Leu_ol at byte {leu_pos}')

sep = b'$$$$'
search_start = max(0, leu_pos - 2000)
preceding = raw.rfind(sep, search_start, leu_pos)
print(f'Preceding $$$$ at byte {preceding}')

if preceding >= 0:
    chunk = raw[preceding:preceding+80]
    print(f'Bytes from $$$$ forward:')
    for i in range(0, len(chunk), 16):
        row = chunk[i:i+16]
        h = ' '.join(f'{b:02x}' for b in row)
        a = ''.join(chr(b) if 32 <= b < 127 else '.' for b in row)
        print(f'  {preceding+i:06x}  {h:<48}  {a}')

# Now try loading just the section from $$$$ to next $$$$
next_sep = raw.find(sep, leu_pos)
print(f'\nNext $$$$ (end of Leu_ol) at byte {next_sep}')

# Extract the Leu_ol record as it appears in the file
record_bytes = raw[preceding + len(sep):next_sep]
record_bytes_with_sep = record_bytes + sep + b'\n'
print(f'Record bytes start: {record_bytes[:40]!r}')

suppl = Chem.SDMolSupplier()
suppl.SetData(record_bytes_with_sep, removeHs=False, sanitize=False)
mols = list(suppl)
print(f'Parse from file bytes: {mols}')
if mols and mols[0]:
    print(f'm_abbr: {mols[0].GetProp("m_abbr")}')

# Count how many records have double CRLF at start
text = raw.decode('utf-8')
records = text.split('$$$$')
double_crlf = [i for i, r in enumerate(records) if r.startswith('\r\n\r\n')]
print(f'\nRecords with double CRLF at start: {len(double_crlf)} — indices: {double_crlf[:10]}')
