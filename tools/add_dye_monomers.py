#!/usr/bin/env python3
"""
Add fluorescent dye maleimide monomers to the pyPept SDF library.

Dyes added:
  5FM      - Fluorescein-5-maleimide  (green, 494/519 nm)
  TAMRA5   - 5-TAMRA maleimide        (orange, 541/568 nm)
  Cy5Mal   - Cyanine5 maleimide       (far-red, 649/670 nm)
  AZDye488 - AZDye 488 maleimide      (green, 499/523 nm)
  AZDye550 - AZDye 550 maleimide      (orange, 555/565 nm)
  AZDye647 - AZDye 647 maleimide      (far-red, 651/672 nm)

All dyes attach via maleimide chemistry (R4 = maleimide_c) to Cys R4 (thiol_s).
"""

import sys
import json
import pathlib
from urllib import request as _req

_REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / 'src'))

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor

SDF_PATH = _REPO / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
rdDepictor.SetPreferCoordGen(True)

MAL_SMARTS = Chem.MolFromSmarts('O=C1C=CC(=O)N1')


def _fetch(url):
    try:
        with _req.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
            p = data['PropertyTable']['Properties'][0]
            return p.get('ConnectivitySMILES') or p.get('CanonicalSMILES') or p.get('IsomericSMILES')
    except Exception as e:
        print(f"    fetch error: {e}")
        return None


def fetch_by_cid(cid):
    return _fetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON")


def fetch_by_name(name):
    return _fetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name.replace(' ', '%20')}/property/CanonicalSMILES/JSON")


def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    return max(frags, key=lambda m: m.GetNumAtoms()) if frags else mol


def add_r4_to_maleimide(mol):
    """Add [4*] substituent to one olefinic C of the maleimide ring."""
    matches = mol.GetSubstructMatches(MAL_SMARTS)
    if not matches:
        return None

    # Pattern: O=C1C=CC(=O)N1
    # indices: 0=O, 1=C(carbonyl), 2=C(olefinic), 3=C(olefinic), 4=C(carbonyl), 5=N
    # Pick match[3] — the olefinic C alpha to the second carbonyl
    target_idx = matches[0][3]

    em = Chem.RWMol(mol)
    r4_idx = em.AddAtom(Chem.Atom(0))          # dummy atom, atomic num 0
    em.GetAtomWithIdx(r4_idx).SetIsotope(4)    # isotope 4 → R4
    em.AddBond(target_idx, r4_idx, Chem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(em)
    except Exception:
        pass

    return em.GetMol()


def load_existing_abbrs():
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    return {mol.GetPropsAsDict().get('m_abbr', '') for mol in suppl if mol}


def register_dye(abbr, display_name, smiles, m_type='linker', m_subtype='dye'):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  FAIL  {abbr}: invalid SMILES")
        return False

    mol = largest_fragment(mol)
    mol = add_r4_to_maleimide(mol)
    if mol is None:
        print(f"  SKIP  {abbr}: no maleimide ring found")
        return False

    try:
        rdDepictor.Compute2DCoords(mol)
    except Exception as e:
        print(f"  WARN  {abbr}: 2D coords failed ({e})")

    mol.SetProp('m_abbr',       abbr)
    mol.SetProp('symbol',       abbr)
    mol.SetProp('m_name',       display_name)
    mol.SetProp('m_type',       m_type)
    mol.SetProp('m_subtype',    m_subtype)
    mol.SetProp('m_Rgroups',    'None,None,None,[H],None,None')
    mol.SetProp('m_chem_types', '4:maleimide_c')

    with open(SDF_PATH, 'a') as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)
        writer.write(mol)
        writer.close()

    print(f"  OK    {abbr} — {display_name} ({mol.GetNumAtoms()} atoms)")
    return True


# (abbr, display_name, fetch_method, identifier)
DYES = [
    ('5FM',      'Fluorescein-5-maleimide',   'cid',  122038),
    ('TAMRA5',   '5-TAMRA maleimide',          'name', 'tetramethylrhodamine-5-maleimide'),
    ('Cy5Mal',   'Cyanine5 maleimide',         'cid',  91757870),
    # AZDye 488 ≈ Alexa Fluor 488 C5 maleimide (sulforhodamine-based, green)
    ('AZDye488', 'AZDye 488 maleimide',        'name', 'alexa fluor 488 maleimide'),
    # AZDye 550 ≈ Cy3 maleimide (cyanine-based, orange/550 nm)
    ('AZDye550', 'AZDye 550 maleimide',        'name', 'cyanine3 maleimide'),
    # AZDye 647 ≈ Alexa Fluor 647 maleimide (cyanine-based, far-red)
    ('AZDye647', 'AZDye 647 maleimide',        'name', 'alexa fluor 647 maleimide'),
]

if __name__ == '__main__':
    existing = load_existing_abbrs()
    print(f"SDF currently has {len(existing)} monomers\n")

    added = 0
    for abbr, display_name, method, ident in DYES:
        if abbr in existing:
            print(f"  SKIP  {abbr}: already in SDF")
            continue

        print(f"Fetching {abbr} ({display_name})…")
        smiles = fetch_by_cid(ident) if method == 'cid' else fetch_by_name(ident)

        if smiles is None:
            print(f"  FAIL  {abbr}: PubChem lookup failed for {method}={ident}")
            continue

        print(f"  SMILES ({len(smiles)} chars): {smiles[:100]}{'…' if len(smiles) > 100 else ''}")
        if register_dye(abbr, display_name, smiles):
            added += 1

    print(f"\nDone — added {added}/{len(DYES)} dyes.")
