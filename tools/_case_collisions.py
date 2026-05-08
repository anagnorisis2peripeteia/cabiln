import pandas as pd, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from rdkit import Chem

CSV = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.csv')
df = pd.read_csv(CSV).set_index('token')

pairs = [('boc', 'Boc'), ('hSer', 'Hser'), ('otbu', 'OtBu')]
for a, b in pairs:
    print(f'--- {a} vs {b} ---')
    for tok in [a, b]:
        if tok in df.index:
            row = df.loc[tok]
            name = row['name']
            smi  = row['input']
            chk  = row['chuckles']
            print(f'  {tok:8s}  name={name}')
            print(f'           input={smi}')
            print(f'           chuckles={chk}')
    print()
