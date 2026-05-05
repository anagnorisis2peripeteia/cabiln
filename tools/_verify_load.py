#!/usr/bin/env python3
"""Verify that the SDF now loads all 1035 records including Leu_ol."""
import pathlib
from rdkit.Chem import PandasTools

SDF = str(pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf')
df = PandasTools.LoadSDF(SDF, removeHs=False, sanitize=False)
print(f'Loaded {len(df)} records')
if 'm_abbr' in df.columns:
    print(f'Leu_ol in m_abbr: {"Leu_ol" in df["m_abbr"].values}')
    leu = df[df['m_abbr'] == 'Leu_ol']
    print(f'Leu_ol rows: {len(leu)}')
    if len(leu) > 0:
        print(f'm_chem_types: {leu.iloc[0].get("m_chem_types", "MISSING")}')
