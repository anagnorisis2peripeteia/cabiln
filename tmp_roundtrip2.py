#!/usr/bin/env python
"""Extended round-trip tests: exotic monomers, functional groups, staples, bicyclics.

Each test has an expected_outcome:
  PASS  - algorithm should handle correctly
  CAP   - terminal cap not in library; check inner residues match
  XLINK - requires crosslink detection (not yet implemented)
"""
import sys, importlib.util
sys.path.insert(0, 'C:/Users/Cameron/repos/pyPept/src')

spec = importlib.util.spec_from_file_location('lr', 'C:/Users/Cameron/repos/pyPept/tools/live_renderer.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

from rdkit import Chem
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence


def cabiln_to_smiles(cabiln):
    seq = Sequence(cabiln)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt='ROMol')
    if romol is None:
        raise ValueError('Assembly returned None')
    return Chem.MolToSmiles(romol)


def cyclic_rotations(seq):
    n = len(seq)
    return [seq[i:] + seq[:i] for i in range(n)]


def cabiln_equiv(a, b):
    def parse(s):
        cyclic = s.startswith('!1-') and s.endswith('-!1')
        inner = s.removeprefix('!1-').removesuffix('-!1') if cyclic else s
        return cyclic, inner.split('-')
    ca, ra = parse(a)
    cb, rb = parse(b)
    if ca != cb or len(ra) != len(rb):
        return False
    return rb in cyclic_rotations(ra) if ca else ra == rb


tests = [
    # ── Bioorthogonal reactive handles ─────────────────────────────────────────
    ('azide AzAla linear',              'G-AzAla-G-L',               'PASS'),
    ('azide AzK linear',                'G-AzK-G-L',                 'PASS'),
    ('azide AzOrn linear',              'G-AzOrn-G-L',               'PASS'),
    ('alkyne Pra linear',               'G-Pra-G-L',                 'PASS'),
    ('allyl AllGly linear',             'G-Gly_allyl-G-L',           'PASS'),
    ('vinyl VinGly linear',             'G-VinGly-G-L',              'PASS'),
    ('azide AzAla cyclic',              '!1-G-AzAla-L-K-!1',        'PASS'),
    ('alkyne Pra cyclic',               '!1-G-Pra-A-F-!1',          'PASS'),
    # ── NHS esters ─────────────────────────────────────────────────────────────
    ('NHSAla linear',                   'G-NHSAla-G-L',             'PASS'),
    ('NHSLys linear',                   'G-NHSLys-G-L',             'PASS'),
    ('NHSDab linear',                   'G-NHSDab-G-L',             'PASS'),
    ('NHSOrn linear',                   'G-NHSOrn-G-L',             'PASS'),
    ('NHSAla cyclic',                   '!1-G-NHSAla-A-L-!1',       'PASS'),
    ('NHSLys cyclic',                   '!1-A-NHSLys-G-P-!1',       'PASS'),
    # ── Maleimides ─────────────────────────────────────────────────────────────
    ('MalAla linear',                   'G-MalAla-G-L',             'PASS'),
    ('MalLys linear',                   'G-MalLys-G-L',             'PASS'),
    ('MalDab linear',                   'G-MalDab-G-L',             'PASS'),
    ('MalOrn linear',                   'G-MalOrn-G-L',             'PASS'),
    ('MalAla cyclic',                   '!1-G-MalAla-A-L-!1',       'PASS'),
    ('MalLys cyclic',                   '!1-A-MalLys-G-P-!1',       'PASS'),
    # ── RCM staple residues pre-crosslink (just alpha-methyl amino acids) ──────
    ('R5 linear pre-staple',            'G-R5-A-L',                 'PASS'),
    ('S5 linear pre-staple',            'G-S5-A-L',                 'PASS'),
    ('R8 linear',                        'G-R8-A-L',                 'PASS'),
    ('S8 linear',                        'G-S8-A-L',                 'PASS'),
    ('R3 linear',                        'G-R3-A-L',                 'PASS'),
    ('S3 linear',                        'G-S3-A-L',                 'PASS'),
    ('R5+S5 two-residue pre-staple',    'G-R5-A-L-A-S5-G',          'PASS'),
    # ── Stapled peptide post-RCM (hydrocarbon crosslink ring) ─────────────────
    # Slot 4 = terminal_alkene; after RCM, slot4-slot4 cyclic alkene bridge
    ('RCM staple i,i+4  S5+S5',        'ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am',   'XLINK'),
    ('RCM staple i,i+7  S5+R8',        'ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am','XLINK'),
    # ── TBMB bicyclic (three Cys each bonded to TBMB crosslinker) ──────────────
    # TBMB not in backbone library; crosslink detection needed for bicyclics
    ('TBMB bicyclic (Cys x3)',
     'ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3', 'XLINK'),
    # ── fmoc / am terminal caps ────────────────────────────────────────────────
    ('fmoc-A-G-L-am',                   'fmoc-A-G-L-am',            'CAP'),
    ('fmoc-K-G-A-F-am',                 'fmoc-K-G-A-F-am',          'CAP'),
    ('fmoc-NHSLys-G-L-am',             'fmoc-NHSLys-G-L-am',       'CAP'),
    ('fmoc-MalAla-G-G-am',             'fmoc-MalAla-G-G-am',       'CAP'),
    ('fmoc-R5-A-S5-G-am pre-staple',   'fmoc-R5-A-S5-G-am',        'CAP'),
    ('ac-A-G-L-am (acetyl cap)',         'ac-A-G-L-am',              'CAP'),
]

ok = fail = xlink = cap_ok = cap_fail = skip = 0
results = []

for label, cabiln, expected in tests:
    # Build SMILES from CABILN
    try:
        smiles = cabiln_to_smiles(cabiln)
    except Exception as e:
        skip += 1
        results.append(('SKIP', label, cabiln, str(e)[:80], ''))
        print(f'SKIP  {label}: CABILN->SMILES: {str(e)[:60]}')
        continue

    if expected == 'XLINK':
        try:
            result, details = mod.smiles_to_cabiln_core(smiles)
            xlink += 1
            results.append(('XLINK', label, cabiln, result, smiles[:60]))
            print(f'XLINK {label}:')
            print(f'      expected: {cabiln}')
            print(f'      got:      {result}')
        except Exception as e:
            xlink += 1
            results.append(('XLINK-ERR', label, cabiln, str(e)[:80], ''))
            print(f'XLINK {label}: algorithm error: {str(e)[:80]}')
        continue

    if expected == 'CAP':
        try:
            result, details = mod.smiles_to_cabiln_core(smiles)
            CAPS = {'fmoc','am','boc','ac','Boc','Fmoc','Ac'}
            def inner_abbrs(s):
                s = s.removeprefix('!1-').removesuffix('-!1')
                return [p for p in s.split('-') if p and p not in CAPS]
            ei = inner_abbrs(cabiln)
            gi = inner_abbrs(result)
            if ei == gi:
                cap_ok += 1
                results.append(('CAP-OK', label, cabiln, result, ''))
                print(f'CAP-OK  {label}: {ei}')
            else:
                cap_fail += 1
                results.append(('CAP-FAIL', label, cabiln, result, smiles[:60]))
                print(f'CAP-FAIL {label}:')
                print(f'         expected inner: {ei}')
                print(f'         got inner:      {gi}')
        except Exception as e:
            cap_fail += 1
            results.append(('CAP-ERR', label, cabiln, str(e)[:80], smiles[:60]))
            print(f'CAP-ERR  {label}: {str(e)[:80]}')
        continue

    # PASS case
    try:
        result, details = mod.smiles_to_cabiln_core(smiles)
    except Exception as e:
        fail += 1
        results.append(('FAIL', label, cabiln, str(e)[:80], smiles[:60]))
        print(f'FAIL  {label}: {str(e)[:80]}')
        continue

    if cabiln_equiv(cabiln, result):
        ok += 1
        print(f'OK    {label}: {cabiln}')
    else:
        fail += 1
        results.append(('FAIL', label, cabiln, result, smiles[:60]))
        print(f'FAIL  {label}:')
        print(f'      expected: {cabiln}')
        print(f'      got:      {result}')

print()
print(f'PASS {ok}/{ok+fail}  |  CAP inner-OK {cap_ok}/{cap_ok+cap_fail}  |  XLINK {xlink} documented  |  SKIP {skip}')

fails = [r for r in results if r[0] in ('FAIL','CAP-FAIL','CAP-ERR')]
if fails:
    print('\nFailures:')
    for r in fails:
        print(f'  [{r[0]}] {r[1]}: got={r[3][:60]}')
