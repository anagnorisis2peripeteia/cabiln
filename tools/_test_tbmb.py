#!/usr/bin/env python3
"""Quick test: can TBMB assemble with three Cys crosslinks?"""

import sys
sys.path.insert(0, "src")

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

# Test 1: TBMB on a pendant chain with three crosslinks
try:
    cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
    print(f"Testing: {cabiln}")
    seq = Sequence(cabiln)
    print(f"  Monomers: {len(seq.s_monomers)}")
    print(f"  Bonds: {seq.s_bonds}")
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt="ROMol")
    if romol:
        from rdkit.Chem import MolToSmiles
        print(f"  SMILES: {MolToSmiles(romol)}")
        print(f"  Atoms: {romol.GetNumAtoms()}")
        print("  OK")
    else:
        print("  FAILED: no molecule")
except Exception as e:
    print(f"  ERROR: {e}")

# Test 2: simple TBMB with just 2 crosslinks (partial use)
try:
    cabiln2 = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2"
    print(f"\nTesting partial: {cabiln2}")
    seq2 = Sequence(cabiln2)
    mol2 = Molecule(seq2)
    romol2 = mol2.get_molecule(fmt="ROMol")
    if romol2:
        from rdkit.Chem import MolToSmiles
        print(f"  SMILES: {MolToSmiles(romol2)}")
        print("  OK")
    else:
        print("  FAILED: no molecule")
except Exception as e:
    print(f"  ERROR: {e}")
