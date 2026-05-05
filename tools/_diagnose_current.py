#!/usr/bin/env python3
"""Diagnose current state of Hsl/Leu_ol in SDF."""
import pathlib, re
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')
print(f'Total records from split: {len(records)}')

# Find Hsl and Leu_ol
hsl_idx = leu_idx = None
for i, rec in enumerate(records):
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m:
        abbr = m.group(1).strip()
        if abbr == 'Hsl':
            hsl_idx = i
        elif abbr == 'Leu_ol':
            leu_idx = i

print(f'Hsl at record index: {hsl_idx}')
print(f'Leu_ol at record index: {leu_idx}')

# Show Hsl's m_subtype / m_chem_types section
if hsl_idx is not None:
    rec = records[hsl_idx]
    m_end = rec.find('M  END')
    props = rec[m_end:m_end+400]
    # Find m_subtype area
    idx = props.find('m_subtype')
    if idx >= 0:
        snippet = props[idx-5:idx+80]
        print(f'\nHsl m_subtype area: {snippet!r}')
    # Find m_chem_types
    idx2 = props.find('m_chem_types')
    if idx2 >= 0:
        snippet2 = props[max(0,idx2-20):idx2+80]
        print(f'Hsl m_chem_types area: {snippet2!r}')
    # Show line endings
    has_crlf = '\r\n' in rec
    has_lf_only = '\n' in rec.replace('\r\n', '')
    print(f'Hsl record: has_crlf={has_crlf}, has_lf_only={has_lf_only}')

# Show Leu_ol's start bytes
if leu_idx is not None:
    rec = records[leu_idx]
    print(f'\nLeu_ol record start (first 60 chars): {rec[:60]!r}')
    has_crlf = '\r\n' in rec
    has_lf_only = '\n' in rec.replace('\r\n', '')
    print(f'Leu_ol record: has_crlf={has_crlf}, has_lf_only={has_lf_only}')

# Find byte offset of Leu_ol record in file
if leu_idx is not None:
    prefix = '$$$$'.join(records[:leu_idx]) + '$$$$'
    offset = len(prefix.encode('utf-8'))
    print(f'\nByte offset where Leu_ol record starts in file: {offset}')
    file_bytes_at_leu = raw[offset:offset+30]
    print(f'File bytes at Leu_ol start: {file_bytes_at_leu!r}')
    # The bytes after $$$$ should be \r\n\r\n     RDKit
    # Find the actual $$$$ that precedes Leu_ol
    sep_bytes = b'$$$$'
    sep_pos = raw.rfind(sep_bytes, 0, offset + 10)
    print(f'Last $$$$ before Leu_ol at byte: {sep_pos}')
    after_sep = raw[sep_pos+4:sep_pos+30]
    print(f'Bytes after that $$$$: {after_sep!r}')

# Test loading
print('\n--- SDMolSupplier test ---')
for sp in [True, False]:
    suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=sp)
    mols = list(suppl)
    total = len(mols)
    leu = sum(1 for m in mols if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Leu_ol')
    hsl = sum(1 for m in mols if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Hsl')
    nones = sum(1 for m in mols if m is None)
    print(f'  strictParsing={sp}: total={total}, nones={nones}, leu={leu}, hsl={hsl}')
