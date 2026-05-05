#!/usr/bin/env python3
"""
Identify exact duplicates in the SDF that should be consolidated via synonyms.
For each duplicate group, pick a canonical entry and make the others synonyms.

Rules for picking canonical:
1. If one token matches a standard 3-letter or 1-letter code, keep that
2. Prefer shorter / more standard name
3. Prefer entry already in CSV (has richer metadata)
"""
import sys, os, csv
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')
CSV_PATH = os.path.join(DATA_DIR, 'monomers.csv')

# Load SDF
print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers")

# Load CSV to see which tokens already have synonyms
csv_rows = {}
with open(CSV_PATH, newline='') as f:
    for row in csv.DictReader(f):
        csv_rows[row['token']] = row

# Group by exact SMILES
groups = defaultdict(list)
for sym in df.index:
    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue
    smi = Chem.MolToSmiles(mol)
    groups[smi].append(sym)

# Find duplicates and assign canonical
print("\n" + "=" * 80)
print("EXACT DUPLICATE GROUPS — PROPOSED CANONICALIZATION")
print("=" * 80)

# Standard AA single-letter codes (these should always be canonical)
STANDARD_AAS = set('GAVLIFPWMSTCYHDNEQKR')

# Preference: in CSV > shorter name > alphabetically first
def pick_canonical(syms):
    # Standard 1-letter AA always wins
    for s in syms:
        if s in STANDARD_AAS:
            return s
    # 3-letter standard (DAla, DVal, etc.) over d-letter codes
    for s in syms:
        if s.startswith('D') and s[1:] in ['Ala','Val','Leu','Ile','Phe','Trp',
                                            'Ser','Thr','Cys','Tyr','His','Asp',
                                            'Glu','Asn','Gln','Lys','Arg']:
            return s
    # Prefer entries in CSV
    in_csv = [s for s in syms if s in csv_rows]
    if in_csv:
        return min(in_csv, key=lambda s: (len(s), s))
    return min(syms, key=lambda s: (len(s), s))

# Collect all synonym assignments
synonym_updates = {}  # canonical -> list of aliases to add
redundant_entries = []  # entries to potentially remove from SDF

total_dupes = 0
for smi, syms in sorted(groups.items(), key=lambda x: -len(x[1])):
    if len(syms) <= 1:
        continue
    total_dupes += 1

    canonical = pick_canonical(syms)
    aliases = [s for s in syms if s != canonical]

    # Check if synonyms already cover these
    existing_syns = csv_rows.get(canonical, {}).get('synonyms', '')
    existing_syn_set = set(s.strip() for s in existing_syns.split(',') if s.strip())

    new_aliases = [a for a in aliases if a not in existing_syn_set]

    if new_aliases:
        synonym_updates[canonical] = {
            'add': new_aliases,
            'existing': existing_syns,
        }
        print(f"\n  Canonical: {canonical}")
        print(f"    Already has synonyms: {existing_syns or '(none)'}")
        print(f"    Add aliases: {', '.join(new_aliases)}")
        print(f"    SMILES: {smi[:70]}")
    else:
        # Already covered
        pass

    redundant_entries.extend(aliases)

print(f"\n\n{'=' * 80}")
print("SUMMARY")
print("=" * 80)
print(f"  Total duplicate groups: {total_dupes}")
print(f"  Groups needing synonym updates: {len(synonym_updates)}")
print(f"  Total redundant SDF entries: {len(redundant_entries)}")
print(f"\n  Redundant entries that could be removed from SDF:")
for r in sorted(redundant_entries):
    in_csv = "  (in CSV)" if r in csv_rows else ""
    print(f"    {r}{in_csv}")

# ═══════════════════════════════════════════════════════════════
# Apply synonym updates to CSV
# ═══════════════════════════════════════════════════════════════
print(f"\n{'=' * 80}")
print("APPLYING SYNONYM UPDATES TO CSV")
print("=" * 80)

rows = []
fieldnames = None
with open(CSV_PATH, newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

updates = 0
for row in rows:
    tok = row['token']
    if tok in synonym_updates:
        existing = row.get('synonyms', '').strip()
        new_aliases = synonym_updates[tok]['add']
        if existing:
            combined = existing + ',' + ','.join(new_aliases)
        else:
            combined = ','.join(new_aliases)
        row['synonyms'] = combined
        updates += 1
        print(f"  {tok}: synonyms = {combined}")

with open(CSV_PATH, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"\n  Updated {updates} CSV entries with new synonyms")
