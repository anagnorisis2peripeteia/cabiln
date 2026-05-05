#!/usr/bin/env python3
"""Check atom counts for affected test sequences."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

def _natoms(biln):
    seq = Sequence(biln)
    mol = Molecule(seq)
    return mol.get_molecule(fmt='ROMol').GetNumAtoms()

sequences = [
    ('C20FA-gGlu-am', 'test_c20fa_gGlu_branch_linear', 32),
    ('C20FA-gGlu-AEEA-am', 'test_c20fa_gGlu_aeea_linear', 42),
    ('ac-K.!1(4,2)-G-am%C20FA-gGlu-AEEA-!1', 'test_lys_sidechain_amide_bond', 58),
]

for cabiln, test_name, old_expected in sequences:
    try:
        actual = _natoms(cabiln)
        delta = actual - old_expected
        print(f"  {test_name}:")
        print(f"    CABILN: {cabiln}")
        print(f"    old expected: {old_expected}, actual: {actual}, delta: {delta:+d}")
    except Exception as e:
        print(f"  {test_name}: ERROR: {e}")

# Retatrutide
print("\n  test_retatrutide:")
retatrutide = (
    "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.!1(4,2)"
    "-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
    "%C20FA-gGlu-AEEA-!1"
)
try:
    actual = _natoms(retatrutide)
    print(f"    old expected: 325, actual: {actual}, delta: {actual - 325:+d}")
except Exception as e:
    print(f"    ERROR: {e}")
