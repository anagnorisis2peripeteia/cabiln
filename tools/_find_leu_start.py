#!/usr/bin/env python3
"""Find the start of Leu_ol's mol block by locating the preceding $$$$."""
import pathlib

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()

# First Leu_ol is in its symbol property; mol block comes earlier.
# Find the mol block start by locating the $$$$ immediately before it.
leu_prop = raw.find(b'Leu_ol')  # in symbol property
# Go back far enough to find the preceding $$$$
search_start = max(0, leu_prop - 2000)
preceding_end = raw.rfind(b'$$$$', search_start, leu_prop)
print(f'Preceding $$$$ at byte {preceding_end}')
print(f'Bytes around it:')
chunk = raw[preceding_end-10:preceding_end+120]
hex_str = ' '.join(f'{b:02x}' for b in chunk)
ascii_str = ''.join(chr(b) if 32 <= b < 127 else '.' for b in chunk)
for i in range(0, len(chunk), 16):
    row = chunk[i:i+16]
    h = ' '.join(f'{b:02x}' for b in row)
    a = ''.join(chr(b) if 32 <= b < 127 else '.' for b in row)
    print(f'  {preceding_end-10+i:06x}  {h:<48}  {a}')
