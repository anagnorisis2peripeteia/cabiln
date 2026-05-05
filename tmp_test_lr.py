#!/usr/bin/env python
"""Test smiles_to_cabiln_core in live_renderer context."""
import sys, importlib.util
sys.path.insert(0, 'C:/Users/Cameron/repos/pyPept/src')
spec = importlib.util.spec_from_file_location('lr', 'C:/Users/Cameron/repos/pyPept/tools/live_renderer.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

tests = [
    ('DKP Gly-LAla', 'C[C@@H]1NC(=O)CNC1=O', '!1-A-G-!1'),  # L-Ala from CABILN round-trip
    ('DKP Gly-DAla', 'O=C1CNC(=O)[C@@H](C)N1', '!1-DAla-G-!1'),  # D-Ala (CIP=R)
    ('cyclo(GGGG)', 'O=C1CNC(=O)CNC(=O)CNC(=O)CN1', '!1-G-G-G-G-!1'),
    ('Arg-containing CP00667', 'C[C@@H](O)[C@@H]1NC(=O)CNC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCCCN)NC1=O', None),
    ('CP00388', 'CC(C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)CNC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC1=O', None),
]

for name, smi, expected in tests:
    try:
        cabiln, details = mod.smiles_to_cabiln_core(smi)
        if expected is None or cabiln == expected:
            status = 'OK'
        else:
            status = f'WRONG (expected {expected})'
        print(f'{name}: {cabiln} [{status}]')
    except Exception as e:
        print(f'{name}: ERROR {e}')
