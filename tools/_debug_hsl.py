#!/usr/bin/env python3
import sys, re, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
from rdkit import Chem
from pyPept.interfaces.monomer_pipeline import (
    pre_activate, _BB_LACTONE_PAT, _BB_COOH_PAT, _BB_ALCOHOL_PAT, _BB_N_PAT
)

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = SDF.read_bytes().decode('latin-1')
records = [r.strip() for r in text.split('$$$$') if r.strip()]

def get_prop(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

for rec in records:
    if get_prop(rec, 'm_abbr') != 'Hsl':
        continue

    rg    = get_prop(rec, 'm_Rgroups')
    ct    = get_prop(rec, 'm_chem_types')
    print(f'Stored rgroups:    {rg}')
    print(f'Stored chem_types: {ct}')

    # Try to get the molblock (everything before first '> ')
    parts = rec.split('> <')
    molblock = parts[0].strip()
    print(f'\nMolblock lines:')
    for l in molblock.split('\n')[:6]:
        print(f'  {repr(l)}')

    # Try parsing with strict=False
    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False, strictParsing=False)
    print(f'\nMolFromMolBlock (strict=False): {mol}')
    if mol:
        smi = Chem.MolToSmiles(mol)
        print(f'SMILES: {smi}')

    # Try via SDMolSupplier with sanitize=False
    suppl = Chem.SDMolSupplier(str(SDF), removeHs=False, sanitize=False)
    for m in suppl:
        if m is None:
            continue
        try:
            p = m.GetPropsAsDict()
        except Exception as e:
            continue
        if p.get('m_abbr', '') == 'Hsl':
            print(f'\nFrom SDMolSupplier:')
            try:
                smi2 = Chem.MolToSmiles(m)
                print(f'SMILES: {smi2}')
                # Test lactone pattern
                mh = Chem.AddHs(m)
                Chem.SanitizeMol(mh)
                print(f'N matches:       {mh.GetSubstructMatches(_BB_N_PAT)}')
                print(f'COOH matches:    {mh.GetSubstructMatches(_BB_COOH_PAT)}')
                print(f'Alcohol matches: {mh.GetSubstructMatches(_BB_ALCOHOL_PAT)}')
                print(f'Lactone matches: {mh.GetSubstructMatches(_BB_LACTONE_PAT)}')
            except Exception as e:
                print(f'Error: {e}')
            break
    break
