#!/usr/bin/env python3
"""
Test smiles_to_cabiln_core on Wikipedia/PubChem-sourced GLP-1 agonist SMILES.
Reports which ones round-trip successfully, prints CABILN for test fixture generation.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import Descriptors

# PubChem / Wikipedia canonical SMILES for each GLP-1 agonist
GLP1_SMILES = {
    'semaglutide': (
        # Wikipedia — MW 4113.6 (C187H291N45O59)
        'CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCC(=O)O)NC(=O)'
        '[C@H](CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@H](NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)'
        'NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(=O)O)'
        'NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)'
        'NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)'
        'NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)NC(=O)[C@@H](N)Cc1c[nH]cn1)'
        '[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)'
        'C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)'
        'N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)C(C)C'
    ),
    'liraglutide': (
        # Wikipedia — MW 3751.3 (C172H265N43O51)
        'CCCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)NCCCC[C@H](NC(=O)[C@H](C)NC(=O)[C@H](C)'
        'NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(C)C)'
        'NC(=O)[C@H](CC1=CC=C(O)C=C1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@@H]'
        '(NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC1=CC=CC=C1)'
        'NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](C)NC(=O)[C@@H](N)'
        'CC1=CN=CN1)[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](CCC(O)=O)C(=O)'
        'N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C)C(=O)'
        'N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)'
        'N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(O)=O)C(O)=O'
    ),
    'exenatide': (
        # Wikipedia — MW 4186.6 (C184H282N50O60S)
        '[H]/N=C(\\N)/NCCC[C@@H](C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1ccccc1)C(=O)'
        'N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2c[nH]c3c2cccc3)'
        'C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)'
        'NCC(=O)N4CCC[C@H]4C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)'
        'C(=O)N5CCC[C@H]5C(=O)N6CCC[C@H]6C(=O)N7CCC[C@H]7C(=O)N[C@@H](CO)C(=O)N)'
        'NC(=O)[C@H](C(C)C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)'
        'NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CCCCN)'
        'NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)'
        'NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](Cc8ccccc8)NC(=O)[C@H]([C@@H](C)O)'
        'NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](Cc9cnc[nH]9)N'
    ),
    'lixisenatide': (
        # PubChem CID 90472060 — MW 4858.6
        'CC[C@H](C)[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)'
        'C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)'
        'NCC(=O)N3CCC[C@H]3C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)'
        'C(=O)N4CCC[C@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCCN)'
        'C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)'
        'C(=O)N[C@@H](CCCCN)C(=O)N)NC(=O)[C@H](CC6=CC=CC=C6)NC(=O)[C@H](CC(C)C)'
        'NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](C(C)C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)'
        'NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCC(=O)N)'
        'NC(=O)[C@H](CCCCN)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)'
        'NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC7=CC=CC=C7)NC(=O)'
        '[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](CC8=CNC=N8)N'
    ),
    'tirzepatide': (
        # PubChem CID 166567236 — MW 4813.5 (C225H348N48O68)
        'CC[C@H](C)[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H]'
        '(CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@@H](C(=O)O)NC(=O)CCCCCCCCCCCCCCCCCCC(=O)O)'
        'C(=O)N[C@@H](C)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](C(C)C)C(=O)'
        'N[C@@H](CCC(=O)N)C(=O)N[C@@H](CC2=CNC3=CC=CC=C32)C(=O)N[C@@H](CC(C)C)'
        'C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C)C(=O)NCC(=O)NCC(=O)N4CCC[C@H]4'
        'C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)C(=O)N5CCC[C@H]5'
        'C(=O)N6CCC[C@H]6C(=O)N7CCC[C@H]7C(=O)N[C@@H](CO)C(=O)N)NC(=O)[C@H](CCCCN)'
        'NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)C(C)(C)NC(=O)[C@H]([C@@H](C)CC)'
        'NC(=O)[C@H](CO)NC(=O)[C@H](CC8=CC=C(C=C8)O)NC(=O)[C@H](CC(=O)O)NC(=O)'
        '[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC9=CC=CC=C9)NC(=O)'
        '[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)'
        'NC(=O)[C@H](CC1=CC=C(C=C1)O)N'
    ),
    'retatrutide': (
        # PubChem CID 171390338 — MW 4731.4
        'CC[C@H](C)C(C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CC=C(C=C1)O)C(=O)'
        'N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)NCC(=O)'
        'NCC(=O)N2CCC[C@H]2C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)'
        'C(=O)N3CCC[C@H]3C(=O)N4CCC[C@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](CO)C(=O)N)'
        'NC(=O)[C@H](CC6=CC=CC=C6)NC(=O)[C@H](C)NC(=O)C(C)(C)NC(=O)[C@H](CCC(=O)N)'
        'NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)COCCOCCNC(=O)CC[C@@H](C(=O)O)'
        'NC(=O)CCCCCCCCCCCCCCCCCCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(=O)O)'
        'NC(=O)[C@H](CC(C)C)NC(=O)[C@@](C)(CC(C)C)NC(=O)C([C@@H](C)CC)NC(=O)'
        '[C@H](CO)NC(=O)[C@H](CC7=CC=C(C=C7)O)NC(=O)[C@H](CC(=O)O)NC(=O)'
        '[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC8=CC=CC=C8)NC(=O)'
        '[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)N)NC(=O)C(C)(C)'
        'NC(=O)[C@H](CC9=CC=C(C=C9)O)N'
    ),
}

print('=== SMILES validation ===\n')
valid = {}
for name, smi in GLP1_SMILES.items():
    mol = Chem.MolFromSmiles(smi)
    if mol:
        mw = Descriptors.MolWt(mol)
        na = mol.GetNumAtoms()
        print(f'{name}: OK  MW={mw:.1f}  atoms={na}')
        valid[name] = smi
    else:
        print(f'{name}: INVALID SMILES')

print()

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

try:
    from tools.live_renderer import smiles_to_cabiln_core
except Exception as e:
    print(f'Import failed: {e}')
    sys.exit(1)

print('=== smiles_to_cabiln round-trip ===\n')
results = {}
for name, smi in valid.items():
    print(f'--- {name} ---')
    try:
        cabiln, details = smiles_to_cabiln_core(smi)
        unknown = [(i, abbr) for i, (abbr, nm, nt) in enumerate(details) if abbr == '?']
        print(f'  CABILN ({len(details)} res): {cabiln}')
        if unknown:
            print(f'  Unknown positions: {unknown[:5]}')
        else:
            results[name] = (smi, cabiln, len(details))
    except Exception as e:
        import traceback
        print(f'  Error: {e}')
        traceback.print_exc()
    print()

print('=' * 60)
print(f'Clean round-trips: {list(results.keys())}')
