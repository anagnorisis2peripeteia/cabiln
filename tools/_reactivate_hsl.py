#!/usr/bin/env python3
"""Patch Hsl's m_chem_types in-place via text substitution (avoids SDWriter encoding issues)."""
import sys, re, pathlib

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

text = SDF.read_bytes().decode('latin-1')
records = text.split('$$$$')

patched = 0
out_records = []
for rec in records:
    if not rec.strip():
        out_records.append(rec)
        continue

    abbr_m = re.search(r'> +<m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    abbr = abbr_m.group(1).strip() if abbr_m else ''

    if abbr == 'Hsl':
        new_ct = '1:backbone_n,2:lactone_c'
        # Replace existing value if present, otherwise insert before end of record
        old_ct_m = re.search(r'(> +<m_chem_types>.*?\r?\n)(.+?)(\r?\n)', rec)
        if old_ct_m:
            old_val = old_ct_m.group(2)
            rec = rec[:old_ct_m.start(2)] + new_ct + rec[old_ct_m.end(2):]
            print(f'Patched Hsl chem_types: {old_val!r} -> {new_ct!r}')
        else:
            # Insert before the trailing blank line / end of record
            eol = '\r\n' if '\r\n' in rec else '\n'
            insert = f'{eol}>  <m_chem_types>  (1) {eol}{new_ct}{eol}'
            rec = rec.rstrip() + insert
            print(f'Inserted Hsl chem_types: {new_ct!r}')
        patched += 1

    out_records.append(rec)

if patched == 0:
    print('ERROR: Hsl not found in SDF')
    sys.exit(1)

SDF.write_bytes('$$$$'.join(out_records).encode('latin-1'))
print('SDF updated.')
