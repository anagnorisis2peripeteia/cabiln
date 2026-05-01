#!/usr/bin/env python3
"""
Tests for pyPept.interfaces.monomer_pipeline — SMARTS-based pre-activation.

Run: python test_leaving_group_pipeline.py
"""
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.warning')

from pyPept.interfaces.monomer_pipeline import (
    pre_activate,
    validate_leaving_group,
)

ALL_20 = [
    ('G', 'Glycine',        'NCC(=O)O'),
    ('A', 'L-Alanine',      'N[C@@H](C)C(=O)O'),
    ('V', 'L-Valine',       'N[C@@H](C(C)C)C(=O)O'),
    ('L', 'L-Leucine',      'N[C@@H](CC(C)C)C(=O)O'),
    ('I', 'L-Isoleucine',   'N[C@@H]([C@@H](C)CC)C(=O)O'),
    ('P', 'L-Proline',      'OC(=O)[C@@H]1CCCN1'),
    ('F', 'L-Phe',          'N[C@@H](Cc1ccccc1)C(=O)O'),
    ('W', 'L-Trp',          'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O'),
    ('M', 'L-Met',          'N[C@@H](CCSC)C(=O)O'),
    ('S', 'L-Ser',          'N[C@@H](CO)C(=O)O'),
    ('T', 'L-Thr',          'N[C@@H]([C@@H](O)C)C(=O)O'),
    ('C', 'L-Cys',          'N[C@@H](CS)C(=O)O'),
    ('Y', 'L-Tyr',          'N[C@@H](Cc1ccc(O)cc1)C(=O)O'),
    ('H', 'L-His',          'N[C@@H](Cc1cnc[nH]1)C(=O)O'),
    ('D', 'L-Asp',          'N[C@@H](CC(=O)O)C(=O)O'),
    ('E', 'L-Glu',          'N[C@@H](CCC(=O)O)C(=O)O'),
    ('N', 'L-Asn',          'N[C@@H](CC(=O)N)C(=O)O'),
    ('Q', 'L-Gln',          'N[C@@H](CCC(=O)N)C(=O)O'),
    ('K', 'L-Lys',          'N[C@@H](CCCCN)C(=O)O'),
    ('R', 'L-Arg',          'N[C@@H](CCCNC(=N)N)C(=O)O'),
]

EXTRA = [
    ('dA', 'D-Alanine',     'N[C@H](C)C(=O)O'),
    ('Aib','Aib',           'CC(N)(C)C(=O)O'),
    ('Hyp','trans-4-Hyp',   'OC(=O)[C@@H]1C[C@@H](O)CN1'),
]

if __name__ == '__main__':
    import sys
    sys.path.insert(0, 'src')

    all_cases = ALL_20 + EXTRA
    passed = failed = 0

    print(f"{'':5} {'AA':<5} {'Name':<20} {'CHUCKLES':<50} {'Leaving groups'}")
    print('-' * 115)

    for aa, name, smi in all_cases:
        ch, lg, _chem_types, err = pre_activate(smi)
        if err:
            print(f"FAIL  {aa:<5} {name:<20} {err}")
            failed += 1
        else:
            lg_str = '  '.join(f'R{s+1}={v}' for s, v in sorted(lg.items()))
            print(f"OK    {aa:<5} {name:<20} {ch:<50} {lg_str}")
            passed += 1

    print()
    print(f"{passed}/{passed+failed} passed")

    print()
    print('-- Validation tests --')
    mol = Chem.AddHs(Chem.MolFromSmiles('N[C@@H](CS)C(=O)O'))
    s_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 16)

    ok, msg = validate_leaving_group(mol, s_idx, '[Cl]')
    assert not ok, "Should reject [Cl] on Cys S"
    print(f"[Cl] correctly rejected: {msg}")

    ok, msg = validate_leaving_group(mol, s_idx, '[H]')
    assert ok, f"Should accept [H] on Cys S: {msg}"
    print(f"[H]  correctly accepted")

    print()
    if failed:
        sys.exit(1)
    print("All tests passed.")
