#!/usr/bin/env python3
"""Test README examples: bracket↔branch roundtrips and molecular equivalence."""

import sys
sys.path.insert(0, "src")

from pyPept.sequence import Sequence, cabiln_to_branch, cabiln_to_bracket
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles
import warnings
warnings.filterwarnings("ignore")

def smiles_for(cabiln):
    seq = Sequence(cabiln)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt="ROMol")
    return MolToSmiles(romol) if romol else None

def test_pair(label, bracket, branch):
    print(f"\n=== {label} ===")
    print(f"  Bracket: {bracket}")
    print(f"  Branch:  {branch}")

    # Test converter roundtrips
    b2r = cabiln_to_branch(bracket)
    r2b = cabiln_to_bracket(branch)
    print(f"  bracket->branch: {b2r}")
    print(f"  branch->bracket: {r2b}")

    # Test molecular equivalence
    smi_bracket = smiles_for(bracket)
    smi_branch = smiles_for(branch)
    print(f"  SMILES bracket: {smi_bracket}")
    print(f"  SMILES branch:  {smi_branch}")

    if smi_bracket == smi_branch:
        print(f"  MATCH: same molecule")
    else:
        print(f"  MISMATCH: different molecules!")

    # Test roundtrips
    rt1 = cabiln_to_bracket(cabiln_to_branch(bracket))
    if rt1 == bracket:
        print(f"  bracket->branch->bracket: OK")
    else:
        print(f"  bracket->branch->bracket: DIFF ({rt1})")

    rt2 = cabiln_to_branch(cabiln_to_bracket(branch))
    if rt2 == branch:
        print(f"  branch->bracket->branch: OK")
    else:
        print(f"  branch->bracket->branch: DIFF ({rt2})")


# README example 1: Lipid linker (lines 143, 146)
test_pair(
    "Lipid linker (K-gGlu-AEEA-C20FA)",
    "ac-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am",
    "ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1",
)

# README example 2: PEG branch (line 136 branch, new bracket)
test_pair(
    "PEG branch",
    "ac-K.[PEG(4,1).am(2,1)]-G-am",
    "ac-K.!1(4,1)-G-am%!1-PEG-am",
)

# README example 3: Fatty acid section (lines 274, 277)
test_pair(
    "Fatty acid lipidation section",
    "ac-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am",
    "ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1",
)

# Test 4: Simple positional branch roundtrip
test_pair(
    "Simple positional branch",
    "ac-C.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-A-G-am",
    "ac-C.(4,4)-A-G-am%gGlu-AEEA(2,1)-C20FA(2,1)",
)

print("\n\nDone.")
