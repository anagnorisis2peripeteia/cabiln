#!/usr/bin/env python3
"""Compare PandasTools vs direct SDMolSupplier behavior."""
import pathlib
from rdkit import Chem
from rdkit.Chem import PandasTools

SDF = str(pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf')

# Direct SDMolSupplier
suppl = Chem.SDMolSupplier(SDF, removeHs=False, sanitize=False)
total = sum(1 for m in suppl)
suppl2 = Chem.SDMolSupplier(SDF, removeHs=False, sanitize=False)
leu = sum(1 for m in suppl2 if m is not None and m.HasProp('m_abbr') and m.GetProp('m_abbr') == 'Leu_ol')
print(f'Direct SDMolSupplier: total={total}, leu_found={leu}')

# PandasTools LoadSDF
df = PandasTools.LoadSDF(SDF, removeHs=False, sanitize=False)
print(f'PandasTools.LoadSDF: rows={len(df)}')
if 'm_abbr' in df.columns:
    print(f'  Leu_ol present: {"Leu_ol" in df["m_abbr"].values}')
    print(f'  Hsl present: {"Hsl" in df["m_abbr"].values}')

# Check PandasTools source
import inspect
print('\nPandasTools.LoadSDF signature:')
print(inspect.signature(PandasTools.LoadSDF))
