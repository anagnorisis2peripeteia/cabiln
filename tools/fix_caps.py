#!/usr/bin/env python3
"""
Fix ACYL-BAD-O caps (R2 on O instead of carbonyl C) and update leaving
groups for all caps based on cap_reactions.yaml.

Two problems being fixed:
1. 16 caps have [2*]OC(=O)... — R2 on ester O creates N-O bonds (wrong).
   Fix: regenerate mol from corrected SMILES with R2 on carbonyl C.
2. All caps use generic [OH]/[H] leaving groups regardless of reaction type.
   Fix: use reagent-specific LGs from cap_reactions.yaml.
"""

import sys, os, csv, re, copy, yaml
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem, Draw

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')
CSV_PATH = os.path.join(DATA_DIR, 'monomers.csv')
YAML_PATH = os.path.join(DATA_DIR, 'cap_reactions.yaml')

# ═══════════════════════════════════════════════════════════════
# Load cap_reactions.yaml
# ═══════════════════════════════════════════════════════════════
with open(YAML_PATH, encoding='utf-8') as f:
    cap_rx = yaml.safe_load(f)

# Map reagent_lg strings to SMILES bracket form for m_Rgroups
LG_MAP = {'Cl': '[Cl]', 'OH': '[OH]', 'H': '[H]', 'I': '[I]', 'Br': '[Br]'}

# ═══════════════════════════════════════════════════════════════
# Load SDF
# ═══════════════════════════════════════════════════════════════
print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers loaded")

# ═══════════════════════════════════════════════════════════════
# Stage 1: Fix ACYL-BAD-O caps
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 1: Fix ACYL-BAD-O caps ===")

fixes = 0
for sym in df.index:
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue
    smi = Chem.MolToSmiles(mol)
    if not smi.startswith('[2*]OC(=O)'):
        continue

    # Fix: replace [2*]OC(=O) with [2*]C(=O)
    fixed_smi = smi.replace('[2*]OC(=O)', '[2*]C(=O)', 1)
    fixed_mol = Chem.MolFromSmiles(fixed_smi)
    if fixed_mol is None:
        print(f"  ERROR: RDKit rejected fixed SMILES for {sym}: {fixed_smi}")
        continue

    # Set isotope on the dummy atom
    for atom in fixed_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetIsotope(2)

    # Generate 2D coords
    AllChem.Compute2DCoords(fixed_mol)

    # Verify the fix
    verify_smi = Chem.MolToSmiles(fixed_mol)
    assert '[2*]C(=O)' in verify_smi or verify_smi.startswith('O=C([2*])'), \
        f"Fix verification failed for {sym}: {verify_smi}"

    df.at[sym, 'ROMol'] = fixed_mol
    fixes += 1
    print(f"  FIXED {sym}: {smi} -> {verify_smi}")

print(f"\n  Fixed {fixes} ACYL-BAD-O caps")

# ═══════════════════════════════════════════════════════════════
# Stage 2: Update leaving groups based on cap_reactions.yaml
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 2: Update leaving groups ===")

# Parse m_Rgroups to list
def parse_rgroups(rg_str):
    if isinstance(rg_str, list):
        return rg_str
    return [None if v.strip() == 'None' else v.strip() for v in rg_str.split(',')]

def format_rgroups(rg_list):
    return ','.join('None' if v is None else v for v in rg_list)

lg_updates = 0
for sym in df.index:
    if sym not in cap_rx:
        continue

    rx = cap_rx[sym]
    reagent_lg = rx.get('reagent_lg')
    if reagent_lg is None:
        continue

    new_lg = LG_MAP.get(str(reagent_lg))
    if new_lg is None:
        print(f"  WARNING: unknown LG '{reagent_lg}' for {sym}")
        continue

    rg = parse_rgroups(df.loc[sym, 'm_Rgroups'])

    # Determine which slot this cap uses
    # R2 caps: slot index 1 (second position)
    # R1 caps: slot index 0 (first position)
    mol = df.loc[sym, 'ROMol']
    smi = Chem.MolToSmiles(mol) if mol else ''

    if '[2*]' in smi:
        slot = 1  # R2
    elif '[1*]' in smi:
        slot = 0  # R1
    else:
        print(f"  WARNING: can't determine slot for {sym}: {smi}")
        continue

    old_lg = rg[slot]
    if old_lg == new_lg:
        continue

    rg[slot] = new_lg
    df.at[sym, 'm_Rgroups'] = format_rgroups(rg)
    lg_updates += 1
    print(f"  {sym}: R{slot+1} LG {old_lg} -> {new_lg} ({rx['reaction']})")

print(f"\n  Updated {lg_updates} leaving groups")

# ═══════════════════════════════════════════════════════════════
# Stage 3: Remove 'issue' flags from cap_reactions.yaml for fixed caps
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 3: Verify final state ===")

# Print all caps with their final SMILES and LGs
caps_in_yaml = sorted(cap_rx.keys())
for sym in caps_in_yaml:
    if sym not in df.index:
        print(f"  WARNING: {sym} in YAML but not in SDF")
        continue
    mol = df.loc[sym, 'ROMol']
    smi = Chem.MolToSmiles(mol) if mol else 'NO_MOL'
    rg = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    # Find active slot
    if '[2*]' in smi:
        lg = rg[1] if len(rg) > 1 else '?'
        slot = 'R2'
    elif '[1*]' in smi:
        lg = rg[0] if len(rg) > 0 else '?'
        slot = 'R1'
    else:
        lg = '?'
        slot = '?'
    rx_type = cap_rx[sym].get('reaction', '?')
    has_issue = 'issue' in cap_rx[sym]
    flag = ' [HAD ISSUE]' if has_issue else ''
    print(f"  {sym:20s} {slot} LG={lg:5s} {rx_type:15s} {smi[:50]}{flag}")

# ═══════════════════════════════════════════════════════════════
# Stage 4: Save fixed SDF
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 4: Save fixed SDF ===")

# Reset index for PandasTools.WriteSDF
df_save = df.copy()
df_save = df_save.reset_index()
df_save = df_save.rename(columns={'index': 'symbol', 'ROMol': 'm_romol'})

# PandasTools needs the mol column named 'ROMol'
if 'm_romol' in df_save.columns:
    df_save = df_save.rename(columns={'m_romol': 'ROMol'})

# Write
backup_path = SDF_PATH + '.bak'
import shutil
shutil.copy2(SDF_PATH, backup_path)
print(f"  Backup saved to {backup_path}")

PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print(f"  SDF written to {SDF_PATH}")

# ═══════════════════════════════════════════════════════════════
# Stage 5: Update CSV
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 5: Update CSV ===")

# Read CSV
rows = []
fieldnames = None
with open(CSV_PATH, newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

csv_updates = 0
for row in rows:
    tok = row['token']
    if tok not in cap_rx:
        continue

    rx = cap_rx[tok]
    reagent_lg = rx.get('reagent_lg')
    if reagent_lg is None:
        continue
    new_lg = LG_MAP.get(str(reagent_lg))
    if new_lg is None:
        continue

    chuckles = row.get('chuckles', '')

    # Fix ACYL-BAD-O in CSV chuckles too
    if '[2*]OC(=O)' in chuckles:
        row['chuckles'] = chuckles.replace('[2*]OC(=O)', '[2*]C(=O)', 1)
        print(f"  FIXED chuckles for {tok}")
        csv_updates += 1

    # Determine slot
    if '[2*]' in chuckles:
        lg_col = 'r2_leaving'
    elif '[1*]' in chuckles:
        lg_col = 'r1_leaving'
    else:
        continue

    old_lg = row.get(lg_col, '')
    if old_lg != new_lg:
        row[lg_col] = new_lg
        csv_updates += 1
        print(f"  {tok}: {lg_col} {old_lg} -> {new_lg}")

# Write CSV
with open(CSV_PATH, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"\n  Updated {csv_updates} CSV entries")

# ═══════════════════════════════════════════════════════════════
# Stage 6: Clean up cap_reactions.yaml — remove issue flags
# ═══════════════════════════════════════════════════════════════
print("\n=== STAGE 6: Update cap_reactions.yaml ===")

# Remove issue flags from fixed caps
yaml_changes = 0
for sym, rx in cap_rx.items():
    if 'issue' in rx:
        del rx['issue']
        yaml_changes += 1
        print(f"  Removed issue flag from {sym}")

with open(YAML_PATH, 'w', encoding='utf-8') as f:
    # Write with custom formatting to preserve structure
    f.write("# Canonical reaction chemistry for each cap monomer.\n")
    f.write("# reaction: acylation | alkylation | sulfonylation | nucleophilic | reductive_amination\n")
    f.write("# reagent_lg: SMILES of the leaving group on the actual reagent\n")
    f.write("# reagent_note: common reagent name\n")
    f.write("\n")
    yaml.dump(cap_rx, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

print(f"\n  Removed {yaml_changes} issue flags")

print("\n" + "=" * 60)
print("ALL DONE")
print("=" * 60)
