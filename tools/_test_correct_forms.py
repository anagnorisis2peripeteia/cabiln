#!/usr/bin/env python3
"""Verify correct bracket/branch equivalence for lipid linker."""

import sys
sys.path.insert(0, "src")

from pyPept.sequence import Sequence, cabiln_to_branch, cabiln_to_bracket
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles
import warnings
warnings.filterwarnings("ignore")

def smiles_for(label, cabiln):
    try:
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt="ROMol")
        smi = MolToSmiles(romol) if romol else "NO MOL"
        print(f"  OK  {label}: {smi[:80]}")
        return smi
    except Exception as e:
        print(f"  FAIL {label}: {e}")
        return None

print("=== CORRECT bracket (1,2) ===")
s1 = smiles_for("bracket",  "ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am")
s2 = smiles_for("branch !1","ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1")
if s1 and s2:
    print(f"  MATCH: {s1 == s2}")

print("\n=== WRONG bracket (2,1) ===")
smiles_for("bracket", "ac-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am")

print("\n=== Converter roundtrips ===")
branch = "ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1"
bracket_from_branch = cabiln_to_bracket(branch)
print(f"  branch -> bracket: {bracket_from_branch}")
s3 = smiles_for("converted bracket", bracket_from_branch)
if s2 and s3:
    print(f"  MATCH branch vs converted: {s2 == s3}")

branch_back = cabiln_to_branch(bracket_from_branch)
print(f"  bracket -> branch: {branch_back}")

bracket_back = cabiln_to_bracket(branch_back)
print(f"  branch -> bracket: {bracket_back}")
print(f"  Full roundtrip OK: {bracket_back == bracket_from_branch}")

print("\n=== Retatrutide bracket (1,2) ===")
smiles_for("retatrutide bracket",
    "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am")

print("\n=== Retatrutide bracket (2,1) ===")
smiles_for("retatrutide bracket",
    "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am")

print("\n=== PEG branch equivalence ===")
s4 = smiles_for("PEG bracket", "ac-K.[PEG(4,1).am(2,1)]-G-am")
s5 = smiles_for("PEG branch",  "ac-K.!1(4,1)-G-am%!1-PEG-am")
if s4 and s5:
    print(f"  MATCH: {s4 == s5}")
