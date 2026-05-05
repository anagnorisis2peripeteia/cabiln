#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from rdkit.Chem import PandasTools

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

df = PandasTools.LoadSDF(SDF_PATH)
print("Columns:", list(df.columns))
print("Shape:", df.shape)
print("\nFirst 3 rows (non-mol):")
for col in df.columns:
    if col != 'ROMol':
        print(f"  {col}: {list(df[col].head(3))}")
