#!/usr/bin/env python3
"""Check which abbrs are missing from loaded records."""
import pathlib, re
from rdkit.Chem import PandasTools

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = SDF.read_bytes().decode('utf-8')
records = text.split('$$$$')

# Get all expected abbrs
expected = set()
for rec in records:
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m:
        expected.add(m.group(1).strip())

print(f'Expected monomers from text split: {len(expected)}')

# Load with PandasTools
df = PandasTools.LoadSDF(str(SDF), removeHs=False, sanitize=False)
print(f'Loaded by PandasTools: {len(df)}')

if 'm_abbr' in df.columns:
    loaded = set(df['m_abbr'].dropna())
    missing = expected - loaded
    print(f'Missing from loaded: {sorted(missing)}')
    extra = loaded - expected
    if extra:
        print(f'Extra in loaded: {sorted(extra)}')
