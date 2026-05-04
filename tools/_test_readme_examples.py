#!/usr/bin/env python3
"""Test every Sequence() call in the README."""
import sys, os, warnings, re, traceback
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
warnings.filterwarnings('ignore')
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit import Chem

readme = os.path.join(os.path.dirname(__file__), '..', 'README.md')
with open(readme, encoding='utf-8') as f:
    text = f.read()

pattern = r"Sequence\(['\"]([^'\"]+)['\"]\)"
matches = re.findall(pattern, text)

EXPECT_FAIL = {
    'C(1,3)-A-A-A-C(1,3)',  # intentional old-BILN rejection demo
}

seen = set()
unique = []
for m in matches:
    if m not in seen:
        seen.add(m)
        unique.append(m)

ok = 0
fail = 0
xfail = 0
errors = []
for seq_str in unique:
    if 'fmt=' in seq_str:
        continue
    expected_fail = seq_str in EXPECT_FAIL
    try:
        seq = Sequence(seq_str)
        mol = Molecule(seq)
        rdmol = mol.get_molecule(fmt='ROMol')
        smi = Chem.MolToSmiles(rdmol)
        if expected_fail:
            print(f'XPASS {seq_str}  (expected failure but passed)', flush=True)
            errors.append((seq_str, 'expected failure but passed'))
            fail += 1
        else:
            print(f'OK  {seq_str}', flush=True)
            ok += 1
    except SystemExit:
        if expected_fail:
            print(f'XFAIL {seq_str}  (expected)', flush=True)
            xfail += 1
        else:
            print(f'ERR {seq_str}', flush=True)
            print(f'    SystemExit (monomer not in library)', flush=True)
            errors.append((seq_str, 'SystemExit - monomer not in library'))
            fail += 1
    except Exception as e:
        if expected_fail:
            print(f'XFAIL {seq_str}  (expected)', flush=True)
            xfail += 1
        else:
            print(f'ERR {seq_str}', flush=True)
            print(f'    {str(e)[:140]}', flush=True)
            errors.append((seq_str, str(e)[:140]))
            fail += 1

print(f'\n{ok} OK, {xfail} expected failures, {fail} FAILED out of {ok+fail+xfail} unique examples', flush=True)
sys.exit(1 if fail else 0)
