#!/usr/bin/env python3
"""Test synonym resolution in sequence assembly."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit import Chem

def assemble(biln):
    seq = Sequence(biln)
    mol = Molecule(seq)
    return Chem.MolToSmiles(mol.get_molecule(fmt='ROMol'))

# Test synonym lookups — these should all resolve to their canonical entry
test_cases = [
    # (using synonym, using canonical, description)
    ('ac-Gly-am',     'ac-G-am',      '3-letter Gly -> G'),
    ('ac-Ala-am',     'ac-A-am',      '3-letter Ala -> A'),
    ('ac-dA-am',      'ac-DAla-am',   'dA -> DAla'),
    ('ac-Nrg-am',     'ac-R-am',      'Nrg -> R (Arg)'),
    ('ButCap-G-am',   'Bua-G-am',     'ButCap -> Bua'),
    ('BzCap-G-am',    'Bz-G-am',      'BzCap -> Bz'),
    ('ac-G-Thi2Ala-am', 'ac-G-Thi-am', 'Thi2Ala -> Thi'),
    ('ac-Dopa-am',    'ac-Tyr_3OH-am', 'Dopa -> Tyr_3OH'),
    ('ac-hCys-am',    'ac-Hcy-am',    'hCys -> Hcy'),
    ('ac-hSer-am',    'ac-Hse-am',    'hSer -> Hse'),
    ('ac-Phe_4OMe-am','ac-Tyr_Me-am', 'Phe_4OMe -> Tyr_Me'),
]

print("=" * 70)
print("SYNONYM RESOLUTION TESTS")
print("=" * 70)

passed = 0
failed = 0

for syn_biln, canon_biln, desc in test_cases:
    try:
        syn_smi = assemble(syn_biln)
        canon_smi = assemble(canon_biln)
        if syn_smi == canon_smi:
            print(f"  PASS: {desc}")
            print(f"        {syn_biln} == {canon_biln}")
            passed += 1
        else:
            print(f"  FAIL: {desc} — different SMILES!")
            print(f"        {syn_biln}: {syn_smi}")
            print(f"        {canon_biln}: {canon_smi}")
            failed += 1
    except Exception as e:
        print(f"  FAIL: {desc} — {e}")
        failed += 1

print(f"\n{'=' * 70}")
print(f"Results: {passed} passed, {failed} failed")
print(f"{'=' * 70}")
