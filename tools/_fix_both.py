#!/usr/bin/env python3
"""Fix both Hsl (missing blank line) and Leu_ol (unknown issue) SDF formatting.

Hsl: the _reactivate_hsl.py script added m_chem_types without the blank line separator
     after m_subtype. Fix: add the missing blank line.

Leu_ol: investigate and fix whatever causes 'Problems encountered parsing data fields'
        with SDMolSupplier strictParsing=True.
"""
import pathlib, re
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

# First, test PandasTools with strictParsing=False to understand Leu_ol
from rdkit.Chem import PandasTools
print('Testing with different SDMolSupplier settings:')
for sp in [True, False]:
    suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=sp)
    found = sum(1 for mol in suppl
                if mol is not None and mol.HasProp('m_abbr') and mol.GetProp('m_abbr') == 'Leu_ol')
    total = sum(1 for mol in suppl)
    # Re-run to count properly
    suppl2 = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False, strictParsing=sp)
    mols = list(suppl2)
    leu_found = sum(1 for m in mols if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Leu_ol')
    total2 = len(mols)
    print(f'  strictParsing={sp}: total={total2}, leu_found={leu_found}')

# Fix Hsl: add missing blank line after m_subtype
fixed_count = 0
for i, rec in enumerate(records):
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m and m.group(1).strip() == 'Hsl':
        # Check if m_subtype is missing blank line before m_chem_types
        if '>  <m_subtype>' in rec:
            # Find m_subtype value section
            sub_m = re.search(r'(>  <m_subtype>.*?\r?\n)(.+?)(\r?\n)(>  <m_chem_types>)', rec, re.DOTALL)
            if sub_m:
                # Check if there's a blank line between m_subtype value and m_chem_types
                between = sub_m.group(3)
                # between should be just \r\n (blank line) but currently it's just \r\n and then >
                # The issue is the value ends with \r\n and then > starts immediately
                orig = sub_m.group(0)
                print(f'Hsl m_subtype section: {orig!r}')
                # Fix: add blank line after value
                eol = '\r\n' if '\r\n' in rec else '\n'
                fixed = orig[:sub_m.end(3) - len(between)] + eol + eol + sub_m.group(4)
                # Actually just: insert blank line before >  <m_chem_types>
                # Current: "natural\r\n>  <m_chem_types>"
                # Target:  "natural\r\n\r\n>  <m_chem_types>"
                old_pattern = sub_m.group(2) + '\r\n' + '>  <m_chem_types>'
                new_pattern = sub_m.group(2) + '\r\n\r\n' + '>  <m_chem_types>'
                if old_pattern in rec:
                    records[i] = rec.replace(old_pattern, new_pattern, 1)
                    print(f'Fixed Hsl: added blank line before m_chem_types')
                    fixed_count += 1
                else:
                    # Try with LF only
                    old_lf = sub_m.group(2) + '\n' + '>  <m_chem_types>'
                    new_lf = sub_m.group(2) + '\n\n' + '>  <m_chem_types>'
                    if old_lf in rec:
                        records[i] = rec.replace(old_lf, new_lf, 1)
                        print(f'Fixed Hsl (LF): added blank line before m_chem_types')
                        fixed_count += 1
                    else:
                        print(f'Could not fix Hsl - pattern not found')
        break

print(f'Fixed {fixed_count} records')
SDF.write_bytes('$$$$'.join(records).encode('utf-8'))
print('Written back.')

# Verify
raw2 = SDF.read_bytes()
text2 = raw2.decode('utf-8')
records2 = text2.split('$$$$')
for rec in records2:
    m = re.search(r'>  <m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m and m.group(1).strip() == 'Hsl':
        m_end = rec.find('M  END')
        props = rec[m_end:]
        print(f'Hsl props after fix: {props!r}')
        break
