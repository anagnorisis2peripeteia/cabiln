#!/usr/bin/env python3
"""Inspect Leu_ol record bytes in detail."""
import pathlib

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()

print(f"First 10 bytes: {raw[:10].hex()} {raw[:10]!r}")

text = raw.decode('utf-8')
records = text.split('$$$$')

leu_rec = None
for i, rec in enumerate(records):
    if 'Leu_ol' in rec and '>  <m_abbr>' in rec:
        leu_rec = rec
        print(f"Leu_ol is record index {i}")
        break

if leu_rec:
    m_end_pos = leu_rec.find('M  END')
    print(f"M  END at char {m_end_pos}")
    props = leu_rec[m_end_pos:]
    print("Properties (repr):")
    print(repr(props))
    print()
    print("Checking for non-ASCII:")
    for j, ch in enumerate(leu_rec):
        if ord(ch) > 127:
            print(f"  char {j}: U+{ord(ch):04X} = {ch!r}")
    print("Done checking.")

    # Also show the mol block header
    header = leu_rec[:200]
    print("\nRecord start (repr):")
    print(repr(header))
