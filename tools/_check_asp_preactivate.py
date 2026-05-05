#!/usr/bin/env python3
"""Check that pre_activate assigns O-attachment carboxyl for standard Asp/Glu."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from rdkit import Chem
from pyPept.interfaces.monomer_pipeline import pre_activate

for name, smi in [("Asp", "N[C@@H](CC(=O)O)C(=O)O"),
                   ("Glu", "N[C@@H](CCC(=O)O)C(=O)O"),
                   ("Asn", "N[C@@H](CC(=O)N)C(=O)O"),
                   ("Ala", "N[C@@H](C)C(=O)O")]:
    result = pre_activate(smi)
    print(f"=== {name} ===")
    print(f"  CHUCKLES: {result.chuckles}")
    print(f"  leaving:  {result.leaving}")
    print(f"  chem_types: {result.chem_types}")

    chuckles = Chem.MolFromSmiles(result.chuckles)
    for a in chuckles.GetAtoms():
        if a.GetAtomicNum() == 0:
            iso = a.GetIsotope()
            nb = list(a.GetNeighbors())[0]
            ct = result.chem_types.get(iso, '?')
            attach = "O-attachment" if nb.GetAtomicNum() == 8 else "C-attachment" if nb.GetAtomicNum() == 6 else "N-attachment"
            print(f"  slot {iso}: dummy -> {nb.GetSymbol()} ({attach}) ct={ct}")
    print()
