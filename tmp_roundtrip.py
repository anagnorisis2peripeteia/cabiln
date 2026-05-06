#!/usr/bin/env python
"""Round-trip tests: CABILN -> SMILES -> CABILN, checking cyclic-rotation equivalence."""
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
        raise ValueError(f'Assembly returned None')
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
    if ca:
        return rb in cyclic_rotations(ra)
    return ra == rb


# ── Test cases ─────────────────────────────────────────────────────────────────
# Each entry: (label, cabiln) or (label, cabiln, synonyms_set)
tests = [
    # ── DKPs: all 20 standard L-amino acids paired with Gly ───────────────────
    ('DKP G-A',                '!1-G-A-!1'),
    ('DKP G-V',                '!1-G-V-!1'),
    ('DKP G-L',                '!1-G-L-!1'),
    ('DKP G-I',                '!1-G-I-!1'),
    ('DKP G-P',                '!1-G-P-!1'),
    ('DKP G-F',                '!1-G-F-!1'),
    ('DKP G-W',                '!1-G-W-!1'),
    ('DKP G-M',                '!1-G-M-!1'),
    ('DKP G-S',                '!1-G-S-!1'),
    ('DKP G-T',                '!1-G-T-!1'),
    ('DKP G-C',                '!1-G-C-!1'),
    ('DKP G-Y',                '!1-G-Y-!1'),
    ('DKP G-N',                '!1-G-N-!1'),
    ('DKP G-Q',                '!1-G-Q-!1'),
    ('DKP G-D',                '!1-G-D-!1'),
    ('DKP G-E',                '!1-G-E-!1'),
    ('DKP G-K',                '!1-G-K-!1'),
    ('DKP G-R',                '!1-G-R-!1'),
    ('DKP G-H',                '!1-G-H-!1'),
    ('DKP F-F',                '!1-F-F-!1'),
    ('DKP P-P',                '!1-P-P-!1'),
    ('DKP W-F',                '!1-W-F-!1'),

    # ── Standard AAs - small/medium cyclic ────────────────────────────────────
    ('cyclo(GGA)',              '!1-G-G-A-!1'),
    ('cyclo(WKP)',              '!1-W-K-P-!1'),
    ('cyclo(GALFVIP)',          '!1-G-A-L-F-V-I-P-!1'),
    ('cyclo(WYHKRDE)',          '!1-W-Y-H-K-R-D-E-!1'),
    ('cyclo(STCMQN)',           '!1-S-T-C-M-Q-N-!1'),
    ('cyclo(KDEKD)',            '!1-K-D-E-K-D-!1'),
    ('cyclo(GALFVIPWYHK)',      '!1-G-A-L-F-V-I-P-W-Y-H-K-!1'),
    ('cyclo(GGGGGGGG)',         '!1-G-G-G-G-G-G-G-G-!1'),
    ('cyclo(PPPP)',             '!1-P-P-P-P-!1'),
    ('cyclo(PPPFFF)',           '!1-P-P-P-F-F-F-!1'),

    # ── Linear peptides ───────────────────────────────────────────────────────
    ('linear GALFVI',           'G-A-L-F-V-I'),
    ('linear WYHKRD',           'W-Y-H-K-R-D'),
    ('linear STCMQN',           'S-T-C-M-Q-N'),
    ('linear GALFVIPWYH',       'G-A-L-F-V-I-P-W-Y-H'),
    ('linear DERKAGPF',         'D-E-R-K-A-G-P-F'),
    ('linear meG-A-L',          'meG-A-L'),

    # ── D-amino acids ─────────────────────────────────────────────────────────
    ('cyclo(DAla-DPhe-G)',      '!1-DAla-DPhe-G-!1'),
    ('cyclo(DLys-DGlu-G-G)',    '!1-DLys-DGlu-G-G-!1'),
    ('cyclo(DAla-A-DPhe-F)',    '!1-DAla-A-DPhe-F-!1'),
    ('cyclo(DVal-G-A)',         '!1-DVal-G-A-!1'),
    ('cyclo(DIle-G-A)',         '!1-DIle-G-A-!1'),
    ('cyclo(DLeu-G-A)',         '!1-DLeu-G-A-!1'),
    ('cyclo(DTrp-G-A)',         '!1-DTrp-G-A-!1'),
    ('cyclo(DTyr-G-A)',         '!1-DTyr-G-A-!1'),
    ('cyclo(DSer-G-A)',         '!1-DSer-G-A-!1'),
    ('cyclo(DThr-G-A)',         '!1-DThr-G-A-!1'),
    ('cyclo(DAsn-G-A)',         '!1-DAsn-G-A-!1'),
    ('cyclo(DGln-G-A)',         '!1-DGln-G-A-!1'),
    ('cyclo(DHis-G-A)',         '!1-DHis-G-A-!1'),
    ('cyclo(DCys-G-A)',         '!1-DCys-G-A-!1'),
    ('cyclo(DPhe-A-G-G)',       '!1-DPhe-A-G-G-!1'),

    # ── N-methylated ──────────────────────────────────────────────────────────
    ('cyclo(meG-A-L-P)',        '!1-meG-A-L-P-!1'),
    ('cyclo(meA-G-G-G)',        '!1-meA-G-G-G-!1'),
    ('cyclo(meF-G-A-L)',        '!1-meF-G-A-L-!1'),
    ('cyclo(meL-meA-G-G)',      '!1-meL-meA-G-G-!1'),
    ('cyclo(meV-G-A-L)',        '!1-meV-G-A-L-!1'),
    ('cyclo(meI-G-A-L)',        '!1-meI-G-A-L-!1'),
    ('cyclo(meK-G-A-L)',        '!1-meK-G-A-L-!1'),
    ('cyclo(meS-G-A-L)',        '!1-meS-G-A-L-!1'),
    ('cyclo(meT-G-A-L)',        '!1-meT-G-A-L-!1'),
    ('cyclo(meW-G-A-L)',        '!1-meW-G-A-L-!1'),
    ('cyclo(meY-G-A-L)',        '!1-meY-G-A-L-!1'),
    ('cyclo(meD-G-A-L)',        '!1-meD-G-A-L-!1'),
    ('cyclo(meE-G-A-L)',        '!1-meE-G-A-L-!1'),

    # ── Non-standard alpha-AAs ────────────────────────────────────────────────
    ('cyclo(Orn-G-A-L)',        '!1-Orn-G-A-L-!1'),
    ('cyclo(Aib-Aib-G-G)',      '!1-Aib-Aib-G-G-!1'),
    ('cyclo(Hyp-G-A-L)',        '!1-Hyp-G-A-L-!1'),
    ('cyclo(Abu-G-A-L)',        '!1-Abu-G-A-L-!1'),
    ('cyclo(Nle-G-A-L)',        '!1-Nle-G-A-L-!1'),
    ('cyclo(Nva-G-A-L)',        '!1-Nva-G-A-L-!1'),
    ('cyclo(Cha-G-A-L)',        '!1-Cha-G-A-L-!1'),
    ('cyclo(Dab-G-A-L)',        '!1-Dab-G-A-L-!1'),
    ('cyclo(Dap-G-A-L)',        '!1-Dap-G-A-L-!1'),
    ('cyclo(Pip-G-A-L)',        '!1-Pip-G-A-L-!1'),
    ('cyclo(Aze-G-A-L)',        '!1-Aze-G-A-L-!1'),
    ('cyclo(Cit-G-A-L)',        '!1-Cit-G-A-L-!1'),
    ('cyclo(Phg-G-A-L)',        '!1-Phg-G-A-L-!1'),
    ('cyclo(Chg-G-A-L)',        '!1-Chg-G-A-L-!1'),
    ('cyclo(Thz-G-A-L)',        '!1-Thz-G-A-L-!1'),
    ('cyclo(Tic-G-A-L)',        '!1-Tic-G-A-L-!1'),
    ('cyclo(Thi-G-A-L)',        '!1-Thi-G-A-L-!1'),
    ('cyclo(Pyr-G-A-L)',        '!1-Pyr-G-A-L-!1'),
    ('cyclo(aIle-G-A-L)',       '!1-aIle-G-A-L-!1'),
    ('cyclo(Hse-G-A-L)',        '!1-Hse-G-A-L-!1'),
    ('cyclo(Hcy-G-A-L)',        '!1-Hcy-G-A-L-!1'),
    ('cyclo(Sec-G-A-L)',        '!1-Sec-G-A-L-!1'),
    ('cyclo(Lys_Ac-G-A-L)',     '!1-Lys_Ac-G-A-L-!1'),

    # ── Modified Phe ──────────────────────────────────────────────────────────
    ('cyclo(Phe_4F-G-G-G)',     '!1-Phe_4F-G-G-G-!1'),
    ('cyclo(Phe_4Cl-A-L-P)',    '!1-Phe_4Cl-A-L-P-!1'),
    ('cyclo(Phe_4OMe-G-G-G)',   '!1-Phe_4OMe-G-G-G-!1', {'Phe_4OMe', 'Tyr_Me'}),
    ('cyclo(Phe_2F-G-A-L)',     '!1-Phe_2F-G-A-L-!1'),
    ('cyclo(Phe_3F-G-A-L)',     '!1-Phe_3F-G-A-L-!1'),
    ('cyclo(Phe_3Cl-G-A-L)',    '!1-Phe_3Cl-G-A-L-!1'),
    ('cyclo(Phe_4Br-G-A-L)',    '!1-Phe_4Br-G-A-L-!1'),
    ('cyclo(Phe_4Me-G-A-L)',    '!1-Phe_4Me-G-A-L-!1'),
    ('cyclo(Phe_4NO2-G-A-L)',   '!1-Phe_4NO2-G-A-L-!1'),
    ('cyclo(Phe_4I-G-A-L)',     '!1-Phe_4I-G-A-L-!1'),
    ('cyclo(Phe_2Me-G-A-L)',    '!1-Phe_2Me-G-A-L-!1'),

    # ── Modified Tyr ──────────────────────────────────────────────────────────
    ('cyclo(Tyr_3OH-G-A-L)',    '!1-Tyr_3OH-G-A-L-!1'),
    ('cyclo(Tyr_3I-G-A-L)',     '!1-Tyr_3I-G-A-L-!1'),
    ('cyclo(Tyr_3NO2-G-A-L)',   '!1-Tyr_3NO2-G-A-L-!1'),

    # ── Naphthylalanines ──────────────────────────────────────────────────────
    ('cyclo(1Nal-G-A-L)',       '!1-1Nal-G-A-L-!1'),
    ('cyclo(2Nal-G-A-L)',       '!1-2Nal-G-A-L-!1'),

    # ── Homo amino acids ──────────────────────────────────────────────────────
    ('cyclo(hPhe-G-A-L)',       '!1-hPhe-G-A-L-!1'),
    ('cyclo(hLys-G-D-E)',       '!1-hLys-G-D-E-!1'),

    # ── Exotic combos ─────────────────────────────────────────────────────────
    ('cyclo(Aib-meG-DPhe-G)',   '!1-Aib-meG-DPhe-G-!1'),
    ('cyclo(Hyp-Orn-meA-G)',    '!1-Hyp-Orn-meA-G-!1'),
    ('cyclo(1Nal-DAla-meF-G)',  '!1-1Nal-DAla-meF-G-!1'),
    ('cyclo(DPhe-A-G-L-V)',     '!1-DPhe-A-G-L-V-!1'),
    ('cyclo(DTrp-meA-G-P)',     '!1-DTrp-meA-G-P-!1'),
    ('cyclo(meF-DLys-G-A-G)',   '!1-meF-DLys-G-A-G-!1'),
    ('cyclo(Chg-meG-DTyr-G)',   '!1-Chg-meG-DTyr-G-!1'),
    ('cyclo(Pip-DAla-meF-G)',   '!1-Pip-DAla-meF-G-!1'),
]


ok = fail = skip = 0
failures = []
skips = []

for entry in tests:
    label, cabiln, *rest = entry
    synonyms: set | None = rest[0] if rest else None
    try:
        smiles = cabiln_to_smiles(cabiln)
    except Exception as e:
        skip += 1
        skips.append((label, cabiln, f'{e}'))
        print(f'SKIP {label}: CABILN->SMILES failed: {e}')
        continue
    try:
        result, details = mod.smiles_to_cabiln_core(smiles)
    except Exception as e:
        fail += 1
        failures.append((label, cabiln, smiles[:60], f'SMILES->CABILN error: {e}'))
        print(f'FAIL {label}: S2C error: {e}')
        continue

    def _accepts(expected, got, synonyms):
        if cabiln_equiv(expected, got):
            return True
        if synonyms:
            for syn in synonyms:
                target = expected.split('-')[1] if expected.startswith('!1-') else expected.split('-')[0]
                alt = got.replace(syn, target)
                if cabiln_equiv(expected, alt):
                    return True
        return False

    if _accepts(cabiln, result, synonyms):
        ok += 1
        note = f' (synonym: {result})' if synonyms and not cabiln_equiv(cabiln, result) else ''
        print(f'OK   {label}: {cabiln}{note}')
    else:
        fail += 1
        failures.append((label, cabiln, smiles[:60], result))
        print(f'FAIL {label}:')
        print(f'     expected: {cabiln}')
        print(f'     got:      {result}')

print(f'\n{ok} OK  {fail} FAIL  {skip} SKIP / {len(tests)} total')
if skips:
    print('\nSkipped:')
    for label, cabiln, reason in skips:
        print(f'  {label}: {reason}')
if failures:
    print('\nFailures detail:')
    for label, cabiln, smiles, got in failures:
        print(f'  {label}: expected={cabiln} got={got}')
