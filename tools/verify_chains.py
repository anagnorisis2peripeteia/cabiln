#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from rdkit import Chem
from rdkit.Chem import PandasTools

df = PandasTools.LoadSDF('src/pyPept/data/monomers.sdf').set_index('symbol')
checks = ['C14FA','C16FA','C18FA','C20FA','C22FA','C14dAcid','C16dAcid','C22dAcid','Myr','Pal','Ste','Ole']
for s in checks:
    mol = df.loc[s, 'ROMol']
    nc = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    print(f"  {s:12s} {nc:2d}C  {Chem.MolToSmiles(mol)[:60]}")
