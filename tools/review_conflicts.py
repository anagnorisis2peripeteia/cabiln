#!/usr/bin/env python3
"""Read-only review: check conflicts in m_chem_types across all SDF monomers."""
import pathlib, re

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

text = SDF.read_bytes().decode('latin-1')
records = [r.strip() for r in text.split('$$$$') if r.strip()]

def get_prop(record, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', record)
    return m.group(1).strip() if m else ''

conflicts = []
missing = []
total = 0

for rec in records:
    total += 1
    abbr  = get_prop(rec, 'm_abbr')
    mtype = get_prop(rec, 'm_type')
    ct    = get_prop(rec, 'm_chem_types')

    if not ct:
        missing.append(abbr or '?')
        continue

    pairs = [p.split(':') for p in ct.split(',') if ':' in p]
    types_set = {v for _, v in pairs}
    has_bb_c = 'backbone_c' in types_set or 'backbone_o' in types_set or 'lactone_c' in types_set
    has_bb_n = 'backbone_n' in types_set or 'cyclooctyne_c' in types_set
    has_backbone = has_bb_n and has_bb_c
    issues = []

    if mtype == 'aa':
        if not has_backbone:
            issues.append('m_type=aa but no complete backbone (N+C/O) detected')
        if has_bb_n and not has_bb_c:
            issues.append('backbone_n without backbone_c or backbone_o')
        if has_bb_c and not has_bb_n:
            issues.append('backbone_c/o without backbone_n')

    if mtype == 'cap' and has_backbone:
        issues.append('m_type=cap but full backbone (N+C) detected -- may be miscategorized')

    # linker: omega-amino acids legitimately have full backbone -- no check

    if len(pairs) > 7:
        issues.append(f'{len(pairs)} R-groups -- unusually many')

    if issues:
        conflicts.append((abbr, mtype, ct, issues))

print(f'Total: {total}  |  Missing chem_types: {len(missing)}  |  Conflicts: {len(conflicts)}')
print()

if missing:
    print('--- STILL MISSING m_chem_types ---')
    print('  ' + ', '.join(missing))
    print()

if conflicts:
    print('--- CONFLICTS ---')
    for abbr, mtype, ct, issues in conflicts:
        print(f'  {abbr} ({mtype}):')
        for iss in issues:
            print(f'    ! {iss}')
        print(f'    chem_types: {ct}')
        print()
