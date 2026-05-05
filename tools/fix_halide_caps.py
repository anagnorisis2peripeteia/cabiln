#!/usr/bin/env python3
"""
Fix 4 hydrocarbon caps by processing their reagent-form SMILES
(with halide LG) through pre_activate.

_Me  → CH3I (methyl iodide)
_Et  → CH3CH2Br (ethyl bromide)
_Bn  → BrCc1ccccc1 (benzyl bromide)
trt  → ClC(c1ccccc1)(c1ccccc1)c1ccccc1 (trityl chloride)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

CAPS = {
    '_Me':  {'reagent': 'CI',                                          'name': 'methyl iodide'},
    '_Et':  {'reagent': 'CCBr',                                       'name': 'ethyl bromide'},
    '_Bn':  {'reagent': 'BrCc1ccccc1',                                'name': 'benzyl bromide'},
    'trt':  {'reagent': 'ClC(c1ccccc1)(c1ccccc1)c1ccccc1',           'name': 'trityl chloride'},
}

print("Testing pre_activate on halide cap reagents:\n")
for sym, info in CAPS.items():
    reagent = info['reagent']
    name = info['name']
    print(f"  {sym} ({name}): {reagent}")
    try:
        result = pre_activate(reagent)
        print(f"    CHUCKLES: {result.chuckles}")
        print(f"    Leaving:  {result.leaving}")
        print(f"    OK")
    except Exception as e:
        print(f"    ERROR: {e}")
    print()

# Now update the SDF
print("=" * 60)
print("Updating SDF entries")
print("=" * 60)

df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')

for sym, info in CAPS.items():
    reagent = info['reagent']

    result = pre_activate(reagent)

    new_mol = Chem.MolFromSmiles(result.chuckles)
    if new_mol is None:
        print(f"  ERROR: can't parse {result.chuckles} for {sym}")
        continue

    AllChem.Compute2DCoords(new_mol)

    max_slot = max(result.leaving.keys()) if result.leaving else 0
    max_slot = max(max_slot, 6)
    rg_parts = [result.leaving.get(i, 'None') for i in range(1, max_slot + 1)]
    new_rgroups = ','.join(rg_parts)

    old_smi = Chem.MolToSmiles(df.loc[sym, 'ROMol'])
    old_rg = df.loc[sym, 'm_Rgroups']

    df.at[sym, 'ROMol'] = new_mol
    df.at[sym, 'm_Rgroups'] = new_rgroups

    print(f"\n  {sym}:")
    print(f"    old CHUCKLES: {old_smi}")
    print(f"    new CHUCKLES: {Chem.MolToSmiles(new_mol)}")
    print(f"    old rgroups:  {old_rg}")
    print(f"    new rgroups:  {new_rgroups}")

# Save
print(f"\nSaving SDF...")
df_save = df.copy().reset_index()
df_save = df_save.rename(columns={'index': 'symbol'})
PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print("  Done")
