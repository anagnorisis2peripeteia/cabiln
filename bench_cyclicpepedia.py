#!/usr/bin/env python3
"""CyclicPepedia SMILES -> CABILN round-trip benchmark.

Downloads CyclicPepedia Peptide_structrure_info.xlsx, extracts SMILES,
runs smiles_to_cabiln_core + cabiln_to_smiles round-trip, reports accuracy
stratified by ring count.

Usage:
    python bench_cyclicpepedia.py [--limit N]
"""
import sys, importlib.util, warnings, json, time, argparse, faulthandler, os
from pathlib import Path

warnings.filterwarnings('ignore')

# Crash forensics: dump C stack to stderr on segfault / fatal signal
faulthandler.enable(file=sys.stderr, all_threads=True)

# ── paths ────────────────────────────────────────────────────────────────────
REPO = Path(__file__).parent.resolve()
sys.path.insert(0, str(REPO / 'src'))

# ── load live_renderer ───────────────────────────────────────────────────────
print('Loading live_renderer...', flush=True)
spec = importlib.util.spec_from_file_location(
    'lr', REPO / 'tools' / 'live_renderer.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
smiles_to_cabiln_core = mod.smiles_to_cabiln_core
print('  OK', flush=True)

from rdkit import Chem
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence


def cabiln_to_smiles(cabiln: str) -> str:
    seq = Sequence(cabiln)
    m = Molecule(seq)
    romol = m.get_molecule(fmt='ROMol')
    if romol is None:
        raise ValueError('Assembly returned None')
    return Chem.MolToSmiles(romol)


def _desalt(mol):
    """Reduce to the parent structure: the single largest organic fragment.

    A '.' in a SMILES means physically separate molecules — counterions
    (acetate, TFA, citrate, mesylate, sulfate, Na/Cl), solvents (chlorobenzene,
    2-methylbutene), waters of crystallisation, scavengers, and metal ions.
    A peptide held together as one entity is covalently linked (disulfide /
    branch) and so is a SINGLE connected component; anything dot-separated is a
    distinct molecule, not part of the peptide being modelled.  The standard
    cheminformatics 'parent' rule therefore applies: keep the largest organic
    component, drop the rest.  Applied symmetrically to source and round-trip,
    so it can never manufacture a false match (within-molecule atom loss stays
    one fragment and is still counted as a mismatch)."""
    frags = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    if len(frags) == 1:
        return mol
    organic = [f for f in frags if any(a.GetAtomicNum() == 6 for a in f.GetAtoms())]
    pool = organic if organic else frags
    return max(pool, key=lambda f: f.GetNumHeavyAtoms())


from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.inchi import MolToInchi
_taut_enum = rdMolStandardize.TautomerEnumerator()


def _normalize(mol):
    """Desalt, iminol-norm, tautomer-canonicalize, then guanidinium-norm last.

    Guanidinium-norm must follow tautomer-canonicalize: the enumerator can flip
    ring-based guanidinium forms (capreomycidine etc.) away from the library
    tautomer.  By running our deterministic norm last it has the final say.
    """
    mol = _desalt(mol)
    mol = mod._s2c_normalize_iminol(mol)
    try:
        mol = _taut_enum.Canonicalize(mol)
    except Exception:
        pass
    mol = mod._s2c_canonicalize_guanidinium(mol)
    return mol


def canonical(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    mol = _normalize(mol)
    return Chem.MolToSmiles(mol)


def canonical_no_stereo(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = _normalize(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def canonical_no_stereo_no_charge(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = _normalize(mol)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        return None


def ring_stratum(n_rings) -> str:
    try:
        n = int(n_rings)
    except (TypeError, ValueError):
        return 'unknown'
    if n <= 1:
        return 'monocyclic (<=1)'
    if n <= 5:
        return 'oligocyclic (2-5)'
    return 'polycyclic (6+)'


# ── download CyclicPepedia xlsx ──────────────────────────────────────────────
import urllib.request

XLSX_URL = ('https://raw.githubusercontent.com/dfwlab/cyclicpepedia'
            '/main/Dataset/Peptide_structrure_info.xlsx')
# Use a local cache next to the script
XLSX_PATH = REPO / 'tools' / 'cyclicpepedia_structure.xlsx'

if not XLSX_PATH.exists():
    print(f'Downloading {XLSX_URL} ...', flush=True)
    try:
        urllib.request.urlretrieve(XLSX_URL, XLSX_PATH)
        print(f'  {XLSX_PATH.stat().st_size // 1024} KB', flush=True)
    except Exception as e:
        print(f'Download failed: {e}', flush=True)
        sys.exit(1)
else:
    print(f'Using cached {XLSX_PATH}', flush=True)

# ── load xlsx ────────────────────────────────────────────────────────────────
try:
    import openpyxl
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'openpyxl', '-q'])
    import openpyxl

wb = openpyxl.load_workbook(XLSX_PATH, read_only=True)
ws = wb.active
rows = list(ws.iter_rows(values_only=True))
headers = list(rows[0])

smiles_col = next(
    (i for i, h in enumerate(headers) if h and str(h).lower() == 'smiles'), None)
rings_col = next(
    (i for i, h in enumerate(headers) if h and 'number_of_rings' in str(h).lower()), None)
if smiles_col is None:
    print('ERROR: no SMILES column found.', flush=True)
    sys.exit(1)
print(f'SMILES col={smiles_col}  rings col={rings_col}', flush=True)

data = rows[1:]

# ── optional row limit ───────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--limit', type=int, default=0)
ap.add_argument('--skip-ids', type=str, default='',
                help='comma-separated peptide IDs to skip (e.g. CP00345,CP01234)')
ap.add_argument('--resume', action='store_true',
                help='resume from cyclicpepedia_bench.partial.json checkpoint')
args = ap.parse_args()
if args.limit:
    data = data[:args.limit]

# crash-loop convergence: auto-skip the row in pending.txt from a prior crash
PENDING_FILE = REPO / 'tools' / 'cyclicpepedia_pending.txt'
SKIP_FILE = REPO / 'tools' / 'cyclicpepedia_skip_ids.txt'
PARTIAL_FILE = REPO / 'tools' / 'cyclicpepedia_bench.partial.json'

skip_ids: set = set()
if args.skip_ids:
    skip_ids.update(s.strip() for s in args.skip_ids.split(',') if s.strip())
if SKIP_FILE.exists():
    skip_ids.update(l.strip() for l in SKIP_FILE.read_text().splitlines() if l.strip())
if PENDING_FILE.exists():
    last_pending = PENDING_FILE.read_text().strip()
    if last_pending and ':' in last_pending:
        idx_s, pid_s = last_pending.split(':', 1)
        print(f'AUTO-SKIP: prior crash at idx={idx_s} id={pid_s}', flush=True)
        skip_ids.add(pid_s.strip())
        with open(SKIP_FILE, 'a') as f:
            f.write(pid_s.strip() + '\n')
    PENDING_FILE.unlink()

if skip_ids:
    print(f'Skipping {len(skip_ids)} known crash-causing IDs', flush=True)

print(f'Rows to process: {len(data)}', flush=True)

# ── run benchmark ────────────────────────────────────────────────────────────
STRATA = ['monocyclic (<=1)', 'oligocyclic (2-5)', 'polycyclic (6+)', 'unknown']
blank = lambda: {'cabiln_ok': 0, 'cabiln_fail': 0,
                 'smiles_exact': 0, 'smiles_no_stereo': 0,
                 'smiles_no_charge': 0, 'inchi_match': 0,
                 'rt_mismatch': 0, 'rt_fail': 0}
totals = blank()
by_stratum: dict = {s: blank() for s in STRATA}
failures: list = []
resume_start = 0

# fold checkpoint back in if --resume
if args.resume and PARTIAL_FILE.exists():
    cp = json.loads(PARTIAL_FILE.read_text())
    totals = cp['summary']
    by_stratum = cp['by_stratum']
    failures = cp.get('failures', [])
    resume_start = cp.get('next_idx', 0)
    print(f'Resumed from checkpoint at idx={resume_start}: '
          f'cabiln_ok={totals["cabiln_ok"]} fail={totals["cabiln_fail"]}', flush=True)

t0 = time.time()

def _checkpoint(next_idx):
    PARTIAL_FILE.write_text(json.dumps({
        'summary': totals, 'by_stratum': by_stratum,
        'failures': failures[:500], 'next_idx': next_idx,
    }, indent=2))

for i, row in enumerate(data):
    if i < resume_start:
        continue
    if i > 0 and i % 500 == 0:
        elapsed = time.time() - t0
        done = totals['cabiln_ok'] + totals['cabiln_fail']
        print(f'  {i}/{len(data)}  cabiln_ok={totals["cabiln_ok"]}  '
              f'fail={totals["cabiln_fail"]}  {elapsed:.1f}s  '
              f'{done/max(elapsed,0.001):.1f}/s', flush=True)
        _checkpoint(i)

    smi = row[smiles_col] if smiles_col < len(row) else None
    if not smi or not isinstance(smi, str) or not smi.strip():
        continue
    smi = smi.strip()
    pid = str(row[0]) if row else f'row_{i}'

    if pid in skip_ids:
        continue

    # forensic marker: written before any C call, deleted after Stage 2
    PENDING_FILE.write_text(f'{i}:{pid}')
    n_rings = row[rings_col] if rings_col is not None and rings_col < len(row) else None
    stratum = ring_stratum(n_rings)
    bucket = by_stratum[stratum]

    # Stage 1: SMILES -> CABILN
    try:
        cabiln, _ = smiles_to_cabiln_core(smi)
        totals['cabiln_ok'] += 1
        bucket['cabiln_ok'] += 1
    except Exception as e:
        totals['cabiln_fail'] += 1
        bucket['cabiln_fail'] += 1
        if len(failures) < 2000:
            failures.append({'id': pid, 'smi': smi, 'stage': 'cabiln',
                             'stratum': stratum, 'err': str(e)[:150]})
        continue

    # Stage 2: CABILN -> SMILES round-trip
    try:
        rt_smi = cabiln_to_smiles(cabiln)
        c_ref, c_rt = canonical(smi), canonical(rt_smi)
        if c_ref == c_rt:
            totals['smiles_exact'] += 1
            bucket['smiles_exact'] += 1
        else:
            ns_ref, ns_rt = canonical_no_stereo(smi), canonical_no_stereo(rt_smi)
            if ns_ref is not None and ns_ref == ns_rt:
                totals['smiles_no_stereo'] += 1
                bucket['smiles_no_stereo'] += 1
                if len(failures) < 2000:
                    failures.append({'id': pid, 'smi': smi, 'cabiln': cabiln,
                                     'rt_smi': rt_smi, 'stage': 'stereo_mismatch',
                                     'stratum': stratum})
            else:
                nc_ref = canonical_no_stereo_no_charge(smi)
                nc_rt = canonical_no_stereo_no_charge(rt_smi)
                if nc_ref is not None and nc_ref == nc_rt:
                    totals['smiles_no_charge'] += 1
                    bucket['smiles_no_charge'] += 1
                    if len(failures) < 2000:
                        failures.append({'id': pid, 'smi': smi, 'cabiln': cabiln,
                                         'rt_smi': rt_smi, 'stage': 'charge_mismatch',
                                         'stratum': stratum})
                else:
                    m_ref = Chem.MolFromSmiles(smi)
                    m_rt = Chem.MolFromSmiles(rt_smi)
                    inchi_ok = False
                    if m_ref and m_rt:
                        try:
                            i_ref = MolToInchi(_desalt(m_ref))
                            i_rt = MolToInchi(_desalt(m_rt))
                            inchi_ok = (i_ref and i_rt and i_ref == i_rt)
                        except Exception:
                            pass
                    if inchi_ok:
                        totals['inchi_match'] += 1
                        bucket['inchi_match'] += 1
                        if len(failures) < 2000:
                            failures.append({'id': pid, 'smi': smi, 'cabiln': cabiln,
                                             'rt_smi': rt_smi, 'stage': 'inchi_match',
                                             'stratum': stratum})
                    else:
                        totals['rt_mismatch'] += 1
                        bucket['rt_mismatch'] += 1
                        if len(failures) < 2000:
                            failures.append({'id': pid, 'smi': smi, 'cabiln': cabiln,
                                             'rt_smi': rt_smi, 'stage': 'rt_mismatch',
                                             'stratum': stratum})
    except Exception as e:
        totals['rt_fail'] += 1
        bucket['rt_fail'] += 1
        if len(failures) < 2000:
            failures.append({'id': pid, 'smi': smi, 'cabiln': cabiln,
                             'stage': 'rt_fail', 'stratum': stratum,
                             'err': str(e)[:150]})

    # row survived — clear forensic marker
    try:
        PENDING_FILE.unlink()
    except FileNotFoundError:
        pass

# ── report ───────────────────────────────────────────────────────────────────
elapsed = time.time() - t0
valid = totals['cabiln_ok'] + totals['cabiln_fail']
pct = lambda n, d: f'{100*n/max(d,1):.1f}%'

print(f'\n{"="*70}')
print(f'CyclicPepedia SMILES -> CABILN -> SMILES Benchmark')
print(f'{"="*70}')
print(f'Total rows: {len(data)}   Valid SMILES: {valid}   Time: {elapsed:.1f}s')
print()
print(f'{"Metric":<28}  {"N":>6}  {"/ valid":>8}  {"/ cabiln_ok":>11}')
print(f'{"-"*58}')
print(f'{"CABILN ok":<28}  {totals["cabiln_ok"]:>6}  {pct(totals["cabiln_ok"],valid):>8}')
print(f'{"CABILN fail":<28}  {totals["cabiln_fail"]:>6}  {pct(totals["cabiln_fail"],valid):>8}')
ok = totals["cabiln_ok"]
cumulative = (totals["smiles_exact"] + totals["smiles_no_stereo"]
              + totals["smiles_no_charge"] + totals["inchi_match"])
print(f'{"RT exact (canonical)":<28}  {totals["smiles_exact"]:>6}  {pct(totals["smiles_exact"],valid):>8}  {pct(totals["smiles_exact"],ok):>11}')
print(f'{"RT match (no stereo)":<28}  {totals["smiles_no_stereo"]:>6}  {pct(totals["smiles_no_stereo"],valid):>8}  {pct(totals["smiles_no_stereo"],ok):>11}')
print(f'{"RT match (no charge)":<28}  {totals["smiles_no_charge"]:>6}  {pct(totals["smiles_no_charge"],valid):>8}  {pct(totals["smiles_no_charge"],ok):>11}')
print(f'{"RT match (InChI)":<28}  {totals["inchi_match"]:>6}  {pct(totals["inchi_match"],valid):>8}  {pct(totals["inchi_match"],ok):>11}')
print(f'{"  => cumulative correct":<28}  {cumulative:>6}  {pct(cumulative,valid):>8}  {pct(cumulative,ok):>11}')
print(f'{"RT graph mismatch":<28}  {totals["rt_mismatch"]:>6}  {"":>8}  {pct(totals["rt_mismatch"],ok):>11}')
print(f'{"RT fail (exception)":<28}  {totals["rt_fail"]:>6}  {"":>8}  {pct(totals["rt_fail"],ok):>11}')

print(f'\n-- By ring count --')
print(f'{"Stratum":<24}  {"N":>5}  {"CABILN%":>8}  {"exact%":>7}  {"no-stereo%":>11}  {"no-charge%":>11}  {"cumul%":>7}')
print(f'{"-" * 80}')
for s in STRATA:
    b = by_stratum[s]
    n = b['cabiln_ok'] + b['cabiln_fail']
    if n == 0:
        continue
    cum = b['smiles_exact'] + b['smiles_no_stereo'] + b['smiles_no_charge'] + b['inchi_match']
    print(f'{s:<24}  {n:>5}  {pct(b["cabiln_ok"],n):>8}  '
          f'{pct(b["smiles_exact"],b["cabiln_ok"]):>7}  '
          f'{pct(b["smiles_no_stereo"],b["cabiln_ok"]):>11}  '
          f'{pct(b["smiles_no_charge"],b["cabiln_ok"]):>11}  '
          f'{pct(cum,b["cabiln_ok"]):>7}')

# failure taxonomy
err_counts: dict = {}
for f in failures:
    k = f.get('err', f.get('stage', '?'))[:80]
    err_counts[k] = err_counts.get(k, 0) + 1
print(f'\n-- Top failure modes (stage + error) --')
stage_counts: dict = {}
for f in failures:
    stage_counts[f.get('stage', '?')] = stage_counts.get(f.get('stage', '?'), 0) + 1
for stage, cnt in sorted(stage_counts.items(), key=lambda x: -x[1]):
    print(f'  {cnt:5d}  [{stage}]')
print()
print(f'-- Top error messages --')
for msg, cnt in sorted(err_counts.items(), key=lambda x: -x[1])[:20]:
    print(f'  {cnt:5d}  {msg}')

out = REPO / 'tools' / 'cyclicpepedia_bench.json'
with open(out, 'w') as f:
    json.dump({'summary': totals, 'by_stratum': by_stratum,
               'total': len(data), 'elapsed': elapsed,
               'failures': failures[:500]}, f, indent=2)
print(f'\nFull results -> {out}')
