#!/usr/bin/env python3
"""
Fix known data errors in the monomer library:
1. C18FA — stored as 16C chain, should be 18C
2. C16dAcid — stored as 14C chain, should be 16C
3. gGlu / D_gGlu — stored identically to E/DGlu, need gamma-linkage
4. mePhe_34diF / D_mePhe_34diF — SMILES is 2F,4Cl which matches neither name
5. mePhe_3Cl / D_mePhe_3Cl — same wrong SMILES as above
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem, Draw

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')
CSV_PATH = os.path.join(DATA_DIR, 'monomers.csv')

# ═══════════════════════════════════════════════════════════════
# Define corrections
# ═══════════════════════════════════════════════════════════════

# Format: symbol -> corrected full SMILES (pre-activation form with all atoms)
# These will be run through the pipeline to generate CHUCKLES

CORRECTIONS = {
    # C18FA: octadecanedioic acid — 18C total (16 CH2 between carbonyls)
    # Was stored as 16C (same as C16FA)
    'C18FA': {
        'chuckles': '[2*]C(=O)CCCCCCCCCCCCCCCCCC(=O)O',
        'reason': 'chain length was 16C, should be 18C (octadecanedioic acid)',
    },
    # C16dAcid: hexadecanedioic acid — 16C total (14 CH2 between carbonyls)
    # Was stored as 14C (same as C14FA)
    'C16dAcid': {
        'chuckles': '[2*]C(=O)CCCCCCCCCCCCCCCC(=O)O',
        'reason': 'chain length was 14C, should be 16C (hexadecanedioic acid)',
    },
    # gGlu: gamma-linked glutamic acid
    # Backbone goes through gamma-COOH (R2), alpha-COOH becomes sidechain (R4)
    # Standard E: [1*]N([3*])[C@@H](CCC([4*])=O)C([2*])=O
    #   -> backbone: N-Calpha-C(=O) (alpha-COOH), sidechain: gamma-COOH
    # gGlu: [1*]N([3*])[C@@H](C([4*])=O)CCC([2*])=O
    #   -> backbone: N-Calpha-CH2-CH2-C(=O) (gamma-COOH), sidechain: alpha-COOH
    'gGlu': {
        'chuckles': '[1*]N([3*])[C@@H](C([4*])=O)CCC([2*])=O',
        'reason': 'was identical to E; gamma-Glu backbone goes through gamma-COOH',
    },
    # D_gGlu: D-gamma-linked glutamic acid
    'D_gGlu': {
        'chuckles': '[1*]N([3*])[C@H](C([4*])=O)CCC([2*])=O',
        'reason': 'was identical to DGlu; D-gamma-Glu backbone goes through gamma-COOH',
    },
    # mePhe_34diF: N-methyl-3,4-difluorophenylalanine
    # Was stored as 2-F,4-Cl (matches neither name)
    'mePhe_34diF': {
        'chuckles': '[1*]N(C)[C@@H](Cc1ccc(F)c(F)c1)C([2*])=O',
        'reason': 'was 2F,4Cl (wrong); corrected to 3,4-difluoro',
    },
    # D_mePhe_34diF: D-form
    'D_mePhe_34diF': {
        'chuckles': '[1*]N(C)[C@H](Cc1ccc(F)c(F)c1)C([2*])=O',
        'reason': 'was 2F,4Cl (wrong); corrected to D-3,4-difluoro',
    },
    # mePhe_3Cl: N-methyl-3-chlorophenylalanine
    # Was stored same as mePhe_34diF (2F,4Cl)
    'mePhe_3Cl': {
        'chuckles': '[1*]N(C)[C@@H](Cc1cccc(Cl)c1)C([2*])=O',
        'reason': 'was identical to mePhe_34diF (2F,4Cl); corrected to 3-chloro',
    },
    # D_mePhe_3Cl: D-form
    'D_mePhe_3Cl': {
        'chuckles': '[1*]N(C)[C@H](Cc1cccc(Cl)c1)C([2*])=O',
        'reason': 'was identical to D_mePhe_34diF; corrected to D-3-chloro',
    },
}

# ═══════════════════════════════════════════════════════════════
# Load and fix SDF
# ═══════════════════════════════════════════════════════════════
print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers loaded\n")

print("=" * 70)
print("FIXING DATA ERRORS IN SDF")
print("=" * 70)

for sym, fix in CORRECTIONS.items():
    if sym not in df.index:
        print(f"  WARNING: {sym} not in SDF, skipping")
        continue

    old_mol = df.loc[sym, 'ROMol']
    old_smi = Chem.MolToSmiles(old_mol) if old_mol else 'NO_MOL'

    new_smi = fix['chuckles']
    new_mol = Chem.MolFromSmiles(new_smi)
    if new_mol is None:
        print(f"  ERROR: RDKit rejected SMILES for {sym}: {new_smi}")
        continue

    # Set isotopes on dummy atoms
    for atom in new_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # Extract isotope from SMILES [n*] pattern
            # The isotope should already be set from parsing [1*], [2*], etc.
            pass

    AllChem.Compute2DCoords(new_mol)
    df.at[sym, 'ROMol'] = new_mol
    print(f"  FIXED {sym}:")
    print(f"    old: {old_smi}")
    print(f"    new: {Chem.MolToSmiles(new_mol)}")
    print(f"    reason: {fix['reason']}")

# ═══════════════════════════════════════════════════════════════
# Save SDF
# ═══════════════════════════════════════════════════════════════
print(f"\nSaving SDF...")
df_save = df.copy().reset_index()
df_save = df_save.rename(columns={'index': 'symbol'})
PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print(f"  SDF written to {SDF_PATH}")

# ═══════════════════════════════════════════════════════════════
# Fix CSV
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("FIXING CSV")
print("=" * 70)

rows = []
fieldnames = None
with open(CSV_PATH, newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

print(f"  CSV columns: {fieldnames}")

csv_fixes = 0
for row in rows:
    tok = row['token']
    if tok not in CORRECTIONS:
        continue

    fix = CORRECTIONS[tok]
    old_chuckles = row.get('chuckles', '')
    row['chuckles'] = fix['chuckles']
    csv_fixes += 1
    print(f"  {tok}: chuckles updated")
    print(f"    old: {old_chuckles}")
    print(f"    new: {fix['chuckles']}")

with open(CSV_PATH, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"\n  Fixed {csv_fixes} CSV entries")

# ═══════════════════════════════════════════════════════════════
# Verify fixes resolved the duplicates
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("VERIFICATION — previously duplicate pairs")
print("=" * 70)

check_pairs = [
    ('C18FA', 'C16FA'),
    ('C16dAcid', 'C14FA'),
    ('gGlu', 'E'),
    ('D_gGlu', 'DGlu'),
    ('mePhe_34diF', 'mePhe_3Cl'),
    ('D_mePhe_34diF', 'D_mePhe_3Cl'),
]

all_resolved = True
for a, b in check_pairs:
    mol_a = df.loc[a, 'ROMol'] if a in df.index else None
    mol_b = df.loc[b, 'ROMol'] if b in df.index else None
    smi_a = Chem.MolToSmiles(mol_a) if mol_a else '?'
    smi_b = Chem.MolToSmiles(mol_b) if mol_b else '?'
    is_dup = (smi_a == smi_b)
    status = "STILL DUPLICATE!" if is_dup else "RESOLVED"
    print(f"  {a} vs {b}: {status}")
    print(f"    {a}: {smi_a}")
    print(f"    {b}: {smi_b}")
    if is_dup:
        all_resolved = False

print(f"\n  {'ALL DUPLICATES RESOLVED' if all_resolved else 'SOME DUPLICATES REMAIN'}")
