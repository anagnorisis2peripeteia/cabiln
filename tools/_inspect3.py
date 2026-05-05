import pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
from rdkit import Chem

sdf = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = sdf.read_bytes().decode('latin-1')
records = [r.strip() for r in text.split('$$$$') if r.strip()]

def get_prop(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

targets = {'deamino_Dab', 'deamino_Lys', 'Hsl'}
for rec in records:
    abbr = get_prop(rec, 'm_abbr')
    if abbr not in targets:
        continue
    mtype = get_prop(rec, 'm_type')
    ct    = get_prop(rec, 'm_chem_types')
    rg    = get_prop(rec, 'm_Rgroups')
    molblock = rec.split('> ')[0].strip()
    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False)
    smi = Chem.MolToSmiles(mol) if mol else '(could not parse)'
    print(f'{abbr} ({mtype}): {smi}')
    print(f'  chem_types: {ct}')
    print(f'  rgroups:    {rg}')
    print()
