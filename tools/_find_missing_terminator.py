#!/usr/bin/env python3
"""Find which record is missing its $$$$ terminator."""
import pathlib, re

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes().decode('utf-8')

# Count $$$$ occurrences
n_terminators = raw.count('$$$$')
print(f'Total $$$$ terminators: {n_terminators}')

records = raw.split('$$$$')
non_empty = [r for r in records if r.strip()]
print(f'Non-empty records: {len(non_empty)}')

def gp(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

# Find records with unusually long M  END to next $$$$ (indicates merged records)
for i, rec in enumerate(records):
    if not rec.strip():
        continue
    # Count M  END lines (should be exactly 1 per record)
    m_end_count = rec.count('M  END')
    if m_end_count != 1:
        abbr = gp(rec, 'm_abbr')
        print(f'Record {i} ({abbr!r}): {m_end_count} M  END lines')
        # Show all abbrs in this merged record
        abbrs = re.findall(r'> +<m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
        print(f'  abbrs found: {abbrs}')

# Also find records where mol block is larger than expected
print('\nRecords with >50 lines before first property:')
for i, rec in enumerate(records):
    if not rec.strip():
        continue
    lines = rec.split('\n')
    # Find first property line
    for j, line in enumerate(lines):
        if line.strip().startswith('>'):
            if j > 40:
                abbr = gp(rec, 'm_abbr')
                print(f'  Record {i} ({abbr}): {j} lines before first property')
            break
