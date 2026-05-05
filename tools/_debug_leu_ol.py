import pathlib, re
from rdkit.Chem import PandasTools, SDMolSupplier

SDF = str(pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf')

# Check if PandasTools loads Leu_ol
df = PandasTools.LoadSDF(SDF, smilesName='SMILES', molColName='Molecule', includeFingerprints=False)
leu = df[df['m_abbr'] == 'Leu_ol'] if 'm_abbr' in df.columns else df[df['symbol'] == 'Leu_ol']
print(f'PandasTools rows total: {len(df)}')
print(f'Leu_ol rows in pandas: {len(leu)}')
if len(leu) == 0:
    print('  -> Leu_ol NOT loaded by PandasTools (mol block parse failure)')

# Compare with SDMolSupplier
suppl = SDMolSupplier(SDF, removeHs=False, sanitize=False)
found = False
for mol in suppl:
    if mol is None:
        continue
    try:
        p = mol.GetPropsAsDict()
        if p.get('m_abbr') == 'Leu_ol':
            found = True
            print(f'\nSDMolSupplier finds Leu_ol: {found}')
            from rdkit import Chem
            try:
                smi = Chem.MolToSmiles(mol)
                print(f'  SMILES: {smi}')
                # Try to sanitize and re-parse
                Chem.SanitizeMol(mol)
                print(f'  Sanitize: OK')
            except Exception as e:
                print(f'  Error: {e}')
            break
    except Exception:
        continue

if not found:
    print('SDMolSupplier also cannot find Leu_ol')
