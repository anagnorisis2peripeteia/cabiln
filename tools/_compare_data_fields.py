#!/usr/bin/env python3
"""Compare data field bytes between Leu_ol (fails) and Hsl (works)."""
import pathlib, re

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

def get_data_fields(rec):
    m_end = rec.find('M  END')
    if m_end < 0:
        return None
    return rec[m_end:]

def hexdump(data, label):
    b = data.encode('utf-8') if isinstance(data, str) else data
    print(f'\n{label} ({len(b)} bytes):')
    for i in range(0, min(len(b), 500), 16):
        row = b[i:i+16]
        h = ' '.join(f'{x:02x}' for x in row)
        a = ''.join(chr(x) if 32 <= x < 127 else '.' for x in row)
        print(f'  {i:04x}  {h:<48}  {a}')

# Find both records
targets = {'Hsl': None, 'Leu_ol': None, 'Lys_ol': None}
for i, rec in enumerate(records):
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m:
        abbr = m.group(1).strip()
        if abbr in targets:
            targets[abbr] = (i, rec)

for name, info in targets.items():
    if info is None:
        print(f'{name}: NOT FOUND')
        continue
    idx, rec = info
    df = get_data_fields(rec)
    hexdump(df, f'{name} (record {idx}) data fields')
    print(f'  repr: {df!r}')
