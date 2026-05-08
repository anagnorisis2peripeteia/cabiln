#!/usr/bin/env python3
"""
Comprehensive monomer library audit:
  1. Exact SMILES duplicates (same canonical free-SMILES, different tokens)
  2. Case-insensitive token collisions
  3. Chirality check — verify CIP at α-carbon matches L-/D- name prefix
     and any explicit (nS,mR) annotations in the name field.

Outputs a report to stdout; exits non-zero if any hard errors are found.
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

CSV_PATH = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data', 'monomers.csv')

# Amino acids (type == 'aa') where L-form has R at α because the side-chain
# has higher CIP priority than COOH.
# Includes: sulfur/selenium directly at β (Cys, Pen, Cys-variants, selenocysteine)
# and β-haloalanines (Cl/Br/F at β outranks O).
# Met is excluded: its sulfur is at γ (CH2-CH2-S), too far to invert α priority.
_L_IS_R_TOKENS = {
    # Cysteine and all its protected/modified variants
    'C', 'DCys', 'dC', 'DC',
    'aCys', 'D_aCys',
    'Cys_All', 'D_Cys_All',
    'Cys_Prg', 'D_Cys_Prg',
    'Cys_tBu', 'D_Cys_tBu',
    'Cys_Trt', 'D_Cys_Trt',
    'Cys_Acm', 'D_Cys_Acm',
    'Cys_Bn', 'D_Cys_Bn',
    'meC', 'D_meC',     # N-methyl Cys
    'aMeCys', 'D_aMeCys',
    # Penicillamine (β,β-dimethyl-Cys)
    'Pen', 'D_Pen',
    # Selenocysteine
    'Sec', 'DSec', 'seC', 'DseC',
    # β-haloalanines (halogen at β beats COOH oxygen in CIP priority)
    'Ala_3F', 'D_Ala_3F',
    'Ala_3Cl', 'D_Ala_3Cl', 'ClAcAla',
    'Ala_3Br', 'D_Ala_3Br',
    'Ala_3I', 'D_Ala_3I',
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def canon(smi):
    """Canonical RDKit SMILES; None if invalid."""
    if not isinstance(smi, str) or not smi.strip():
        return None
    m = Chem.MolFromSmiles(smi.strip())
    return Chem.MolToSmiles(m) if m else None

def strip_stereo(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)

def alpha_cip(smi, token=''):
    """
    Return the CIP code at the putative α-carbon of an amino acid.
    The α-carbon is the first chiral centre that bears N, C=O, and one
    other non-H substituent.  Returns 'S', 'R', or None.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    # Look for: N-C(*)-C(=O) pattern — the α-carbon
    patt = Chem.MolFromSmarts('[NX3][C@@H,C@H,C@@,C@]([CX3](=O))')
    if patt is None:
        return None
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None
    # Take the first match; index 1 is the α-carbon
    alpha_idx = matches[0][1]
    atom = mol.GetAtomWithIdx(alpha_idx)
    return atom.GetPropsAsDict().get('_CIPCode')

def parse_name_stereo(name):
    """
    Extract expected L/D and any (nS,mR,...) from the name string.
    Uses the FIRST L-/D- word-boundary match so "D-4-methyl-L-proline"
    correctly returns 'D' (not 'L' from the embedded "L-proline").
    Returns ('L'|'D'|None, list_of_(position,code) or []).
    """
    if not isinstance(name, str):
        return None, []
    # Find first occurrence of word-boundary L- or D-
    first_l = re.search(r'\bL-', name)
    first_d = re.search(r'\bD-', name)
    ld = None
    if first_l and first_d:
        ld = 'L' if first_l.start() < first_d.start() else 'D'
    elif first_l:
        ld = 'L'
    elif first_d:
        ld = 'D'
    elif re.match(r'l-', name, re.IGNORECASE):
        ld = 'L'
    # Explicit CIP annotations like (2S,3R) or (S) or (R)
    cip_pat = re.findall(r'\((\d*[SR](?:,\d*[SR])*)\)', name, re.IGNORECASE)
    explicit = []
    for grp in cip_pat:
        for part in grp.split(','):
            part = part.strip()
            m = re.match(r'^(\d*)([SR])$', part, re.IGNORECASE)
            if m:
                pos = int(m.group(1)) if m.group(1) else None
                code = m.group(2).upper()
                explicit.append((pos, code))
    return ld, explicit


def expected_alpha_cip(ld, token, has_sulfur_in_sc=None):
    """
    Given 'L' or 'D', return the expected CIP at α for this token.
    Cys/Met types are special: L→R, D→S.
    """
    if ld is None:
        return None
    l_is_r = token in _L_IS_R_TOKENS
    if ld == 'L':
        return 'R' if l_is_r else 'S'
    else:
        return 'S' if l_is_r else 'R'


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print('Loading CSV...')
df = pd.read_csv(CSV_PATH)
print(f'  {len(df)} monomers\n')

# ---------------------------------------------------------------------------
# 1. Exact SMILES duplicates
# ---------------------------------------------------------------------------
print('=' * 60)
print('1. EXACT SMILES DUPLICATES (same canonical free-SMILES)')
print('=' * 60)

# Compute canonical SMILES for each row
df['_canon'] = df['input'].apply(canon)

valid = df[df['_canon'].notna()].copy()
dup_mask = valid.duplicated(subset='_canon', keep=False)
dup_df = valid[dup_mask].sort_values('_canon')

exact_dup_groups = []
for canon_smi, grp in dup_df.groupby('_canon'):
    tokens = grp['token'].tolist()
    names  = grp['name'].tolist()
    types  = grp['type'].tolist()
    exact_dup_groups.append((tokens, names, canon_smi))

if exact_dup_groups:
    print(f'Found {len(exact_dup_groups)} duplicate groups:\n')
    for tokens, names, smi in sorted(exact_dup_groups, key=lambda x: x[0]):
        print(f"  tokens: {tokens}")
        print(f"  names:  {names}")
        print(f"  SMILES: {smi[:80]}")
        print()
else:
    print('  None found.\n')

# ---------------------------------------------------------------------------
# 2. Core duplicates (strip stereo, same flat structure)
# ---------------------------------------------------------------------------
print('=' * 60)
print('2. CORE-ONLY DUPLICATES (same connectivity, different stereo or token)')
print('=' * 60)

valid['_flat'] = valid['_canon'].apply(lambda s: strip_stereo(s) if s else None)
flat_groups = []
for flat, grp in valid[valid['_flat'].notna()].groupby('_flat'):
    if len(grp) < 2:
        continue
    # Only flag if NOT already in exact dup list
    canon_set = set(grp['_canon'].tolist())
    if len(canon_set) == 1:
        continue   # all are exact dups, already reported
    tokens = grp['token'].tolist()
    names  = grp['name'].tolist()
    canons = grp['_canon'].tolist()
    flat_groups.append((tokens, names, canons, flat))

if flat_groups:
    print(f'Found {len(flat_groups)} core-only duplicate groups (differ only in stereo):\n')
    for tokens, names, canons, flat in sorted(flat_groups):
        print(f"  tokens: {tokens}")
        print(f"  names:  {names}")
        for t, c in zip(tokens, canons):
            print(f"    {t:12s}  {c[:80]}")
        print()
else:
    print('  None found.\n')

# ---------------------------------------------------------------------------
# 3. Case-insensitive token collisions
# ---------------------------------------------------------------------------
print('=' * 60)
print('3. CASE-INSENSITIVE TOKEN COLLISIONS')
print('=' * 60)

df['_token_lower'] = df['token'].str.lower()
ci_dup = df[df.duplicated(subset='_token_lower', keep=False)].sort_values('_token_lower')
ci_groups = []
for lower, grp in ci_dup.groupby('_token_lower'):
    rows = grp[['token', 'name', 'input', '_canon']].values.tolist()
    ci_groups.append((lower, rows))

if ci_groups:
    print(f'Found {len(ci_groups)} case-insensitive collision groups:\n')
    for lower, rows in ci_groups:
        print(f"  '{lower}':")
        for tok, name, smi, can in rows:
            same = '(same structure)' if can == rows[0][3] else '*** DIFFERENT STRUCTURE ***'
            print(f"    {tok:12s}  {name:40s}  {same}")
        print()
else:
    print('  None found.\n')

# ---------------------------------------------------------------------------
# 4. Chirality audit for amino acids
# ---------------------------------------------------------------------------
print('=' * 60)
print('4. CHIRALITY AUDIT (amino acids with L-/D- or explicit CIP in name)')
print('=' * 60)

aa_df = df[(df['type'] == 'aa') & df['_canon'].notna()].copy()

errors = []
warnings = []
ok_count = 0
skipped = 0

for _, row in aa_df.iterrows():
    tok  = row['token']
    name = row['name']
    smi  = row['_canon']

    ld, explicit = parse_name_stereo(name)
    actual_cip   = alpha_cip(smi, tok)

    # Check explicit (nS/nR) annotations in name vs actual CIP
    if explicit:
        # Only use position-2 annotation if present, else single-code
        pos2 = [code for pos, code in explicit if pos == 2 or pos is None]
        if pos2:
            exp = pos2[0]
            if actual_cip and actual_cip != exp:
                errors.append({
                    'token': tok, 'name': name, 'expected': exp,
                    'actual': actual_cip, 'reason': 'explicit CIP annotation in name'
                })
            elif actual_cip == exp:
                ok_count += 1
            else:
                skipped += 1
            continue

    # Check L-/D- prefix
    if ld is None:
        skipped += 1
        continue

    exp = expected_alpha_cip(ld, tok)
    if actual_cip is None:
        # No chiral centre found — could be Gly or a flat molecule
        skipped += 1
        continue

    if actual_cip != exp:
        errors.append({
            'token': tok, 'name': name, 'expected': exp,
            'actual': actual_cip, 'reason': f'{ld}-prefix convention'
        })
    else:
        ok_count += 1

print(f'OK:      {ok_count}')
print(f'Skipped: {skipped}  (no L/D prefix and no explicit CIP annotation)')
print(f'ERRORS:  {len(errors)}\n')

if errors:
    print('--- CHIRALITY MISMATCHES ---\n')
    for e in errors:
        print(f"  {e['token']:12s}  {e['name']:45s}  expected={e['expected']}  actual={e['actual']}")
        print(f"               reason: {e['reason']}")
else:
    print('All audited amino acids have correct alpha-carbon CIP.\n')

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print('\n' + '=' * 60)
print('SUMMARY')
print('=' * 60)
print(f'  Total monomers:              {len(df)}')
print(f'  Exact duplicate groups:      {len(exact_dup_groups)}')
print(f'  Core-only duplicate groups:  {len(flat_groups)}')
print(f'  Case-collision groups:       {len(ci_groups)}')
print(f'  Chirality errors:            {len(errors)}')

has_errors = len(errors) > 0 or any(
    rows[0][3] != rows[i][3]
    for _, rows in ci_groups
    for i in range(1, len(rows))
)
sys.exit(1 if has_errors else 0)
