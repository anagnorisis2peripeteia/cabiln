#!/usr/bin/env python3
"""Change m_type=cap -> aa for deamino_Dab and deamino_Lys (pipeline detects full backbone)."""
import sys, re, pathlib

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
TARGETS = {'deamino_Dab', 'deamino_Lys'}

text = SDF.read_bytes().decode('latin-1')
records = text.split('$$$$')

patched = 0
out_records = []
for rec in records:
    if rec.strip():
        abbr_m = re.search(r'> +<m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
        abbr = abbr_m.group(1).strip() if abbr_m else ''
        if abbr in TARGETS:
            mtype_m = re.search(r'(> +<m_type>.*?\r?\n)(.+?)(\r?\n)', rec)
            if mtype_m:
                old = mtype_m.group(2)
                rec = rec[:mtype_m.start(2)] + 'aa' + rec[mtype_m.end(2):]
                print(f'{abbr}: m_type {old!r} -> "aa"')
                patched += 1
    out_records.append(rec)

SDF.write_bytes('$$$$'.join(out_records).encode('latin-1'))
print(f'\nPatched {patched}/{len(TARGETS)} entries.')
