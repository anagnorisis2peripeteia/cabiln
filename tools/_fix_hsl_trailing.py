#!/usr/bin/env python3
"""Fix Hsl's missing blank line after m_chem_types value (before $$$$)."""
import pathlib, re
from rdkit import Chem
from rdkit.Chem import PandasTools

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')
print(f'Records from split: {len(records)}')

# Find Hsl
hsl_idx = None
for i, rec in enumerate(records):
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m and m.group(1).strip() == 'Hsl':
        hsl_idx = i
        break

print(f'Hsl at index: {hsl_idx}')
rec = records[hsl_idx]

# Show the tail of Hsl's record
print(f'Hsl record tail: {rec[-80:]!r}')

# Fix: ensure blank line before the end of the record (before $$$$)
# The record should end with \r\n\r\n (blank line)
# Currently ends with \r\n (single CRLF)
fixed = False
if rec.endswith('\r\n\r\n'):
    print('Hsl already has trailing blank line - no fix needed')
elif rec.endswith('\r\n'):
    records[hsl_idx] = rec + '\r\n'
    print('Fixed: added trailing \\r\\n blank line to Hsl')
    fixed = True
elif rec.endswith('\n\n'):
    print('Hsl already has trailing blank line (LF) - no fix needed')
elif rec.endswith('\n'):
    records[hsl_idx] = rec + '\n'
    print('Fixed: added trailing \\n blank line to Hsl')
    fixed = True
else:
    print(f'Unexpected record ending: {rec[-10:]!r}')

if fixed:
    new_text = '$$$$'.join(records)
    SDF.write_bytes(new_text.encode('utf-8'))
    print(f'Written. New file size: {SDF.stat().st_size} bytes')

# Verify
print('\n--- Verification ---')
for sp in [True, False]:
    suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=sp)
    mols = list(suppl)
    leu = sum(1 for m in mols if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Leu_ol')
    hsl = sum(1 for m in mols if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Hsl')
    print(f'  SDMolSupplier strictParsing={sp}: total={len(mols)}, leu={leu}, hsl={hsl}')

df = PandasTools.LoadSDF(str(SDF), removeHs=False, sanitize=False)
print(f'  PandasTools: total={len(df)}')
if 'm_abbr' in df.columns:
    for t in ['Leu_ol', 'Hsl']:
        rows = df[df['m_abbr'] == t]
        print(f'    {t}: {len(rows)} rows', end='')
        if len(rows) > 0:
            ct = rows.iloc[0].get('m_chem_types', 'MISSING')
            print(f', m_chem_types={ct!r}', end='')
        print()
