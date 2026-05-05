#!/usr/bin/env python3
"""Check how Leu_ol's record starts and test isolated parsing."""
import pathlib
from rdkit import Chem

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
raw = SDF.read_bytes()
text = raw.decode('utf-8')
records = text.split('$$$$')

for i, rec in enumerate(records):
    if 'Leu_ol' in rec and '>  <m_abbr>' in rec:
        print(f'Leu_ol record index: {i}')
        print(f'Record starts with: {rec[:60]!r}')
        print(f'Leading chars: {[hex(ord(c)) for c in rec[:6]]}')

        # Test isolated
        test_sdf = (rec + '$$$$\n').encode('utf-8')
        suppl = Chem.SDMolSupplier()
        suppl.SetData(test_sdf, removeHs=False, sanitize=False)
        mols = list(suppl)
        print(f'Isolated parse: {mols}')
        if mols and mols[0]:
            print(f'm_abbr: {mols[0].GetProp("m_abbr")}')
        break
