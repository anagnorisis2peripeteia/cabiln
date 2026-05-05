import sys
sys.path.insert(0, 'src')
from rdkit import Chem
from pyPept.interfaces.reaction_library import infer_chem_type

SDF_PATH = 'src/pyPept/data/monomers.sdf'
targets = {'D_meD', 'D_meE', 'Ald', 'Lys_For'}
suppl = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
for mol in suppl:
    if mol is None:
        continue
    props = mol.GetPropsAsDict()
    abbr = props.get('m_abbr', '')
    if abbr not in targets:
        continue
    smi = Chem.MolToSmiles(mol)
    rg = props.get('m_Rgroups', '')
    ct = props.get('m_chem_types', '')
    print(f'{abbr}: {smi}')
    print(f'  RG={rg}')
    print(f'  CT={ct}')
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            iso = atom.GetIsotope()
            if iso > 0:
                attach = list(atom.GetNeighbors())[0]
                inferred = infer_chem_type(mol, attach.GetIdx(), iso)
                print(f'  slot{iso}: attach={attach.GetSymbol()}[{attach.GetIdx()}] inferred={inferred}')
    print()
