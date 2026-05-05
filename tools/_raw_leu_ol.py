#!/usr/bin/env python3
"""Dump raw bytes of Leu_ol's property section and surrounding records."""
import pathlib, re

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()

# Find "Leu_ol" as bytes
leu_pos = raw.find(b'Leu_ol')
print(f'First occurrence of b"Leu_ol" at byte {leu_pos}')

# Show 200 bytes before and 400 after
start = max(0, leu_pos - 200)
chunk = raw[start:leu_pos + 400]
print(f'\nRaw bytes ({len(chunk)} bytes from pos {start}):')
print(chunk)
print()

# Show as hex for non-printable / non-ASCII
print('Hex dump:')
for i in range(0, len(chunk), 16):
    row = chunk[i:i+16]
    hex_part = ' '.join(f'{b:02x}' for b in row)
    ascii_part = ''.join(chr(b) if 32 <= b < 127 else '.' for b in row)
    print(f'  {start+i:06x}  {hex_part:<48}  {ascii_part}')
