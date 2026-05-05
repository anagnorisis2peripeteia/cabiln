#!/usr/bin/env python3
import sys, re, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
from rdkit import Chem
from pyPept.interfaces.monomer_pipeline import pre_activate, _BB_LACTONE_PAT, ActivationError

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

def restore_plain(chuckles, rg_str):
    rgroups = [r.strip() for r in rg_str.split(',')]
    s = chuckles
    for i, lg in enumerate(rgroups):
        slot = i + 1
        if lg.lower() not in ('none', ''):
            s = re.sub(r'\[' + str(slot) + r'\*\]', lg, s)
    return s

suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False)
for mol in suppl:
    if mol is None:
        continue
    try:
        p = mol.GetPropsAsDict()
    except Exception as e:
        continue
    if p.get('m_abbr', '') != 'Hsl':
        continue

    rg = p.get('m_Rgroups', '')
    print(f'Found Hsl. rgroups={rg}')

    chuckles = Chem.MolToSmiles(mol)
    print(f'CHUCKLES: {chuckles}')

    plain = restore_plain(chuckles, rg)
    print(f'Plain:    {plain}')

    # Test lactone pattern on plain SMILES mol
    m2 = Chem.MolFromSmiles(plain)
    print(f'MolFromSmiles: {m2}')
    if m2:
        m2h = Chem.AddHs(m2)
        hits = m2h.GetSubstructMatches(_BB_LACTONE_PAT)
        print(f'Lactone hits on plain+H: {hits}')

    # Now try pre_activate
    try:
        result = pre_activate(plain)
        print(f'pre_activate chem_types: {result.chem_types}')
    except ActivationError as e:
        print(f'ActivationError: {e}')
    except Exception as e:
        print(f'Exception: {type(e).__name__}: {e}')

    break
