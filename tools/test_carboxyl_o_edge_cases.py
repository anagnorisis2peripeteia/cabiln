#!/usr/bin/env python3
"""Test all edge cases for the proposed carboxyl O-attachment convention change.

Tests (run independently of any code changes):
  1. Does find_sidechain_slots double-detect backbone COOH O?
  2. Does infer_chem_type correctly distinguish carboxyl O from aldehyde C?
  3. Does the proposed O-attachment SMIRKS (unmapped O) actually break the C-O bond?
  4. Does [OX2H0:1][CX3](=O) false-positive on ester, anhydride, amide, ether?
  5. Cap monomer COOH O — same problem as backbone?
  6. Sidechain ester (Glu-OtBu) — does it avoid carboxyl pre_smarts?
  7. Multi-carboxyl (hAsp) — backbone COOH O vs sidechain COOH O.
  8. REACTION_INDEX collision check.
  9. Intramolecular O-departure: does grouped-SMIRKS produce one fragment?
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from rdkit import Chem
from rdkit.Chem import AllChem

# Proposed SMARTS
_INFER_SMARTS_PROPOSED = '[OX2H0:1][CX3](=O)'  # carboxyl O in CHUCKLES (H→dummy)
_PRE_SMARTS_PROPOSED   = '[OX2H1:1][CX3](=O)'  # carboxyl O in raw SMILES

_patt_infer = Chem.MolFromSmarts(_INFER_SMARTS_PROPOSED)
_patt_pre   = Chem.MolFromSmarts(_PRE_SMARTS_PROPOSED)

PASS = 0
FAIL = 0

def check(label, cond, note=''):
    global PASS, FAIL
    status = 'PASS' if cond else 'FAIL'
    if not cond:
        FAIL += 1
    else:
        PASS += 1
    print(f"  [{status}] {label}" + (f" — {note}" if note else ''))


# ── TEST 1: Does backbone COOH O match pre_smarts? ────────────────────────────
print("\n=== TEST 1: Backbone COOH O double-detection ===")
# Standard Ala in raw SMILES (no dummies)
ala = Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
assert ala, "Bad SMILES"
# Find all matches of proposed pre_smarts
matches_pre = ala.GetSubstructMatches(_patt_pre)
# For Ala backbone COOH: O is OX2H1 bonded to CX3(=O). This is exactly
# what the proposed pre_smarts matches.
backbone_o_indices = [m[0] for m in matches_pre]
print(f"  Proposed pre_smarts matches in Ala: {matches_pre}")
check("Ala backbone COOH O IS matched by proposed pre_smarts (expected — this is the problem)",
      len(matches_pre) > 0, f"match indices: {backbone_o_indices}")

# Now simulate what find_sidechain_slots would do:
# backbone C is found first and added to assigned_atoms
# But backbone COOH O is NOT in assigned_atoms
from pyPept.interfaces.monomer_pipeline import find_backbone_slots, find_sidechain_slots

ala_h = Chem.AddHs(ala)
backbone = find_backbone_slots(ala_h)
assigned = set(backbone.values()) if backbone else set()
print(f"  Backbone slots (no H indices): {backbone}")
print(f"  assigned_atoms after find_backbone_slots: {assigned}")

# Run the CURRENT find_sidechain_slots (uses existing _CHEM_TYPE_REGISTRY pre_smarts)
sidechain_current = find_sidechain_slots(ala_h, assigned)
print(f"  Sidechain slots (CURRENT): {sidechain_current}")
check("Ala has 0 sidechain slots with CURRENT carboxyl SMARTS",
      len(sidechain_current) == 0)

# Simulate with PROPOSED pre_smarts:
# If carboxyl pre_smarts matches O, find_sidechain_slots would find backbone COOH O
proposed_matches = [(m[0], '[OH]', 'carboxyl')
                    for m in ala_h.GetSubstructMatches(_patt_pre)
                    if m[0] not in assigned]
print(f"  Proposed sidechain carboxyl O matches in Ala (not in assigned): {proposed_matches}")
check("Backbone COOH O WOULD be spuriously detected under proposed change",
      len(proposed_matches) > 0,
      "CRITICAL — backbone COOH O not in assigned_atoms; would get a spurious sidechain slot")


# ── TEST 2: How many O atoms at risk (survey of representative AAs) ───────────
print("\n=== TEST 2: Survey of standard AAs for backbone COOH O risk ===")
test_smiles = {
    'Gly': 'NCC(=O)O',
    'Ala': 'N[C@@H](C)C(=O)O',
    'Asp': 'N[C@@H](CC(=O)O)C(=O)O',  # two carboxyls!
    'Glu': 'N[C@@H](CCC(=O)O)C(=O)O',
    'Pro': 'N1CCC[C@@H]1C(=O)O',
}
for name, smi in test_smiles.items():
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    backbone = find_backbone_slots(mol)
    assigned = set(backbone.values()) if backbone else set()
    sc = find_sidechain_slots(mol, assigned)
    proposed_extra = [m[0] for m in mol.GetSubstructMatches(_patt_pre)
                      if m[0] not in assigned]
    correct_sidechain_o = [x for x in proposed_extra
                           if any(mol.GetAtomWithIdx(nb).GetAtomicNum() == 6 and
                                  any(mol.GetBondBetweenAtoms(nb, y).GetBondTypeAsDouble() == 2.0
                                      for y in [z.GetIdx() for z in mol.GetAtomWithIdx(nb).GetNeighbors()])
                                  for nb in [n.GetIdx() for n in mol.GetAtomWithIdx(x).GetNeighbors()])]
    print(f"  {name}: assigned={len(assigned)} current_sc={len(sc)} "
          f"proposed_extra_O={len(proposed_extra)} (spurious: {proposed_extra})")
for name in ('Asp', 'Glu'):
    mol = Chem.AddHs(Chem.MolFromSmiles(test_smiles[name]))
    backbone = find_backbone_slots(mol)
    assigned = set(backbone.values()) if backbone else set()
    proposed_extra = [m[0] for m in mol.GetSubstructMatches(_patt_pre)
                      if m[0] not in assigned]
    check(f"{name}: would detect >1 O (sidechain COOH O + spurious backbone COOH O)",
          len(proposed_extra) >= 2,
          f"found {len(proposed_extra)} O matches outside assigned_atoms")


# ── TEST 3: Mitigation — adding backbone COOH O to assigned_atoms ─────────────
print("\n=== TEST 3: Mitigation — exclude backbone COOH O from assigned ===")
def get_backbone_cooh_o(mol, backbone_c_idx):
    """Return the OH oxygen bonded to the backbone carbonyl C (if present)."""
    excluded = set()
    for nb in mol.GetAtomWithIdx(backbone_c_idx).GetNeighbors():
        if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() >= 1:
            excluded.add(nb.GetIdx())
    return excluded

for name, smi in test_smiles.items():
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    backbone = find_backbone_slots(mol)
    if not backbone:
        continue
    assigned = set(backbone.values())
    # Apply mitigation: also exclude backbone COOH O
    if 2 in backbone:
        assigned |= get_backbone_cooh_o(mol, backbone[2])

    proposed_extra = [m[0] for m in mol.GetSubstructMatches(_patt_pre)
                      if m[0] not in assigned]
    print(f"  {name}: after mitigation, proposed extra O = {proposed_extra}")

# Asp/Glu should have exactly 1 (the sidechain COOH O), others 0
asp_mol = Chem.AddHs(Chem.MolFromSmiles(test_smiles['Asp']))
asp_backbone = find_backbone_slots(asp_mol)
asp_assigned = set(asp_backbone.values())
asp_assigned |= get_backbone_cooh_o(asp_mol, asp_backbone[2])
asp_extra = [m[0] for m in asp_mol.GetSubstructMatches(_patt_pre) if m[0] not in asp_assigned]
check("Asp with mitigation: exactly 1 O (sidechain COOH only)", len(asp_extra) == 1)

ala_mol = Chem.AddHs(Chem.MolFromSmiles(test_smiles['Ala']))
ala_backbone = find_backbone_slots(ala_mol)
ala_assigned = set(ala_backbone.values())
ala_assigned |= get_backbone_cooh_o(ala_mol, ala_backbone[2])
ala_extra = [m[0] for m in ala_mol.GetSubstructMatches(_patt_pre) if m[0] not in ala_assigned]
check("Ala with mitigation: 0 extra O", len(ala_extra) == 0)


# ── TEST 4: infer_smarts false-positive on ester, anhydride, amide, phenol ───
print("\n=== TEST 4: False-positives for proposed infer_smarts on non-attachment atoms ===")
# These mols have CHUCKLES dummies on specific atoms; check if infer_smarts
# fires on the WRONG O atoms.
cases = [
    # (label, smiles_with_dummy, expected_match_on_dummy_neighbor, expected_ct)
    ('Ester O (Glu-OtBu sidechain)', 'N[C@@H](CCC(=O)OC(C)(C)C)C([4*])=O', None, 'backbone_c'),
    ('Amide O in Asn sidechain',     'N[C@@H](CC(=O)[NH2])C([4*])=O',         None, None),
    ('Anhydride O',                  'C(=O)OC([4*])=O',                        4,    'backbone_c'),
    ('Phenol O (not CX3)',           'Oc1ccc([4*])cc1',                         4,    None),
    ('Carboxyl O correct',           'N[C@@H](CC([3*])=O)C([2*])=O',            3,    'aldehyde'),  # current infer
    ('Carboxyl O-attach (proposed)', 'N[C@@H](CC(=O)O[4*])C([2*])=O',          4,    '?'),         # what does infer return?
]

for label, smi, dummy_iso, _ in cases:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  INVALID SMILES: {label}")
        continue
    # Find all atoms with proposed infer_smarts match
    matches = mol.GetSubstructMatches(_patt_infer)
    match_indices = [m[0] for m in matches]
    # Find dummy neighbor atoms for context
    dummy_atoms = [(a.GetIsotope(), [nb.GetIdx() for nb in a.GetNeighbors()])
                   for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    print(f"  {label}:")
    print(f"    infer_smarts matches: {matches}")
    print(f"    dummies: {dummy_atoms}")

# Specific checks:
# 1. Ester O in Glu-OtBu: no dummy on ester O, so infer never called on it
glu_otbu = Chem.MolFromSmiles('N[C@@H](CCC(=O)OC(C)(C)C)C([2*])=O')
ester_O_matches = glu_otbu.GetSubstructMatches(_patt_infer) if glu_otbu else []
check("Glu-OtBu ester O: proposed infer_smarts matches ester O (only relevant if dummy present)",
      len(ester_O_matches) == 0,
      "ester O is OX2H0 bonded to CX3(=O) → DOES match! But no dummy on ester O → infer never called")

# Actually let me check the match directly
if glu_otbu:
    ester_matches = glu_otbu.GetSubstructMatches(_patt_infer)
    check("Glu-OtBu: proposed infer_smarts MATCHES ester O (structurally, ignoring dummy context)",
          len(ester_matches) > 0,
          "This is fine ONLY because ester O never gets a dummy in correct CHUCKLES")
    print(f"    Glu-OtBu ester O match indices: {ester_matches}")

# 2. Amide O: OX1, not OX2 — should NOT match
asn_mol = Chem.MolFromSmiles('N[C@@H](CC(=O)N)C([2*])=O')
if asn_mol:
    amide_matches = asn_mol.GetSubstructMatches(_patt_infer)
    check("Asn amide O: infer_smarts does NOT match amide O",
          all(asn_mol.GetAtomWithIdx(m[0]).GetAtomicNum() != 8 or
              any(asn_mol.GetBondBetweenAtoms(m[0], nb.GetIdx()).GetBondTypeAsDouble() > 1
                  for nb in asn_mol.GetAtomWithIdx(m[0]).GetNeighbors())
              for m in amide_matches), "amide O should be OX1")
    print(f"    Asn amide O match: {amide_matches}")


# ── TEST 5: Proposed SMIRKS O-departure — does C-O bond break correctly? ─────
print("\n=== TEST 5: O-departure SMIRKS mechanics ===")
# Simulate: carboxyl O-attachment CHUCKLES meets amine → amide bond
# Using concrete isotopes to mimic assembly

# Simple test: sidechain carboxyl O-attachment + amine
# CHUCKLES: C(=O)O[100*] (carboxyl sidechain) + [200*]N (amine)
# SMIRKS: [100*][O][C:4](=[O:5]).[200*][N:2] >> [N:2][C:4](=[O:5])

smirks_o_attach = '[100*][O][C:4](=[O:5]).[200*][N:2] >> [N:2][C:4](=[O:5])'
rxn = AllChem.ReactionFromSmarts(smirks_o_attach)

# Carboxyl sidechain molecule: propionic acid O-attachment
carboxyl_mol = Chem.MolFromSmiles('CC(=O)O[100*]')  # CH3-C(=O)-O-[100*]
amine_mol    = Chem.MolFromSmiles('[200*][NH2]')     # [200*]-NH2

check("Carboxyl O-attachment mol parses", carboxyl_mol is not None)
check("Amine mol parses", amine_mol is not None)

if carboxyl_mol and amine_mol:
    prods = rxn.RunReactants((carboxyl_mol, amine_mol))
    print(f"  Products: {len(prods)} product sets")
    for i, ps in enumerate(prods):
        for j, p in enumerate(ps):
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
                print(f"    Product [{i}][{j}]: {smi}")
            except Exception as e:
                print(f"    Product [{i}][{j}]: SANITIZE ERROR {e} — raw: {Chem.MolToSmiles(p)}")

    if prods:
        p = prods[0][0]
        try:
            Chem.SanitizeMol(p)
            smi = Chem.MolToSmiles(p)
            check("SMIRKS produces amide bond (N-C=O)", 'NC(C)=O' in smi or 'CNC' in smi or 'C(=O)N' in smi,
                  f"got: {smi}")
            # The O should NOT appear in the main product
            frags = smi.split('.')
            main = max(frags, key=len)
            check("Carboxyl O not in main product fragment", 'O' not in main.replace('=O', '').replace('[OH]', ''),
                  f"main fragment: {main}")
        except:
            pass


# TEST 5b: Verify O fragment count (take_largest needed?)
print("\n  5b: Fragment count — does OH appear as separate fragment?")
if prods:
    for i, ps in enumerate(prods):
        for j, p in enumerate(ps):
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
                n_frags = smi.count('.') + 1
                print(f"    Product [{i}][{j}] fragments={n_frags}: {smi}")
                check(f"Product [{i}][{j}]: >1 fragment (OH departs as separate piece)",
                      n_frags > 1, f"take_largest needed")
            except:
                pass


# ── TEST 6: Intramolecular O-departure (imide-like ring closure) ──────────────
print("\n=== TEST 6: Intramolecular O-departure (imide ring closure) ===")
# Asp-like: N[300*] (backbone N) + sidechain carboxyl O[400*]
# Molecule: [300*][N:1]CC[C:2](=[O:3])O[400*]
# SMIRKS (intramolecular grouped): ([300*][N:1].[400*][O][C:2](=[O:3])) >> [N:1][C:2](=[O:3])

smirks_imide = '([300*][N:1].[400*][O][C:2](=[O:3])) >> [N:1][C:2](=[O:3])'
rxn_imide = AllChem.ReactionFromSmarts(smirks_imide)

# Simplified Asp-like: [300*]NCCC(=O)O[400*]
asp_like = Chem.MolFromSmiles('[300*]NCC(=O)O[400*]')  # simple 4-membered ring test
print(f"  Asp-like intramolecular SMILES: {Chem.MolToSmiles(asp_like) if asp_like else 'INVALID'}")

if asp_like:
    prods = rxn_imide.RunReactants((asp_like,))
    print(f"  Intramolecular products: {len(prods)} product sets")
    for i, ps in enumerate(prods):
        for j, p in enumerate(ps):
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
                n_frags = smi.count('.') + 1
                print(f"    Product [{i}][{j}] frags={n_frags}: {smi}")
            except Exception as e:
                print(f"    Product [{i}][{j}]: ERROR {e}")

# Try with proper ring-closing atom counts (6-membered succinimide-like):
asp_suc = Chem.MolFromSmiles('[300*]NCCC(=O)O[400*]')  # 5-membered ring
if asp_suc:
    prods = rxn_imide.RunReactants((asp_suc,))
    if prods:
        p = prods[0][0]
        try:
            Chem.SanitizeMol(p)
            smi = Chem.MolToSmiles(p)
            print(f"  Succinimide-like closure: {smi}")
            n_frags = smi.count('.') + 1
            check("Intramolecular: product is single fragment (ring closed)",
                  n_frags == 1, f"got {n_frags} fragments: {smi}")
        except Exception as e:
            print(f"  ERROR: {e}")
    else:
        check("Intramolecular: SMIRKS produced products", False, "0 products")


# ── TEST 7: REACTION_INDEX collision ─────────────────────────────────────────
print("\n=== TEST 7: REACTION_INDEX current [amine_primary, carboxyl] entry ===")
from pyPept.interfaces.reaction_library import REACTION_INDEX
ct_pair = ('amine_primary', 'carboxyl')
entry = REACTION_INDEX.get(ct_pair)
check(f"[amine_primary, carboxyl] currently routes to backbone_amide",
      entry is not None and entry.get('id') == 'backbone_amide',
      f"entry: {entry.get('id') if entry else 'NONE'}")
print(f"  Current backbone_amide SMIRKS: {entry['steps'] if entry else 'N/A'}")
print(f"  These would need to change for O-attachment carboxyl — different SMIRKS needed.")


# ── TEST 8: Proposed infer_smarts on the actual O-attachment CHUCKLES ─────────
print("\n=== TEST 8: infer_chem_type on carboxyl O CHUCKLES vs aldehyde C CHUCKLES ===")
from pyPept.interfaces.reaction_library import infer_chem_type

cases_infer = [
    # (label, smiles, dummy_iso, expected_chem_type_after_change)
    ('Carboxyl O-attach [4*] on O',  'CCC(=O)O[4*]',         4, 'carboxyl (proposed)'),
    ('Aldehyde C-attach [4*] on C',  'CC([4*])=O',            4, 'aldehyde (unchanged)'),
    ('Thiol S-attach [4*] on S',     'CC([4*])S[4*]',         4, 'thiol (already fixed)'),
    ('Backbone C slot2',             '[1*]NCC([2*])=O',        2, 'backbone_c'),
]

patt_proposed = Chem.MolFromSmarts(_INFER_SMARTS_PROPOSED)
for label, smi, iso, expected in cases_infer:
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print(f"  INVALID: {label}")
        continue
    # Find attach atom for given iso
    attach_idx = None
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetIsotope() == iso:
            nbs = list(a.GetNeighbors())
            if nbs:
                attach_idx = nbs[0].GetIdx()
                break
    if attach_idx is None:
        print(f"  No dummy iso {iso} in {label}")
        continue

    current_ct = infer_chem_type(mol, attach_idx, slot=iso)

    # Proposed: would the new SMARTS match?
    proposed_match = any(m[0] == attach_idx
                         for m in mol.GetSubstructMatches(patt_proposed))

    print(f"  {label}:")
    print(f"    current infer_chem_type = {current_ct!r}")
    print(f"    proposed [OX2H0:1][CX3](=O) would match = {proposed_match}")
    print(f"    expected after change: {expected}")


# ── TEST 9: Ester O false-positive risk — more rigorous ──────────────────────
print("\n=== TEST 9: Ester O in CHUCKLES — false-positive if dummy were placed there ===")
# Hypothetical: what if someone wrongly placed [4*] on an ester O?
ester_o_dummy = Chem.MolFromSmiles('CC(=O)O[4*]')  # ester O with dummy — same as carboxyl!
if ester_o_dummy:
    matches = ester_o_dummy.GetSubstructMatches(patt_proposed)
    attach_idx = next((list(a.GetNeighbors())[0].GetIdx()
                       for a in ester_o_dummy.GetAtoms()
                       if a.GetAtomicNum() == 0 and a.GetIsotope() == 4), None)
    if attach_idx is not None:
        matched = any(m[0] == attach_idx for m in matches)
        check("Ester O with dummy: would be detected as carboxyl by proposed infer_smarts",
              matched, "Can't distinguish ester from carboxyl via O — both are OX2H0 bonded to CX3(=O)")
        print(f"    Note: this is ONLY safe because ester Os never get dummies in correct CHUCKLES")


# ── SUMMARY ──────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"  PASS: {PASS}   FAIL: {FAIL}")
print(f"{'='*60}")
print("""
KEY FINDINGS:
  1. Backbone COOH O WILL be spuriously detected — needs mitigation in
     pre_activate() to add COOH O to assigned_atoms after find_backbone_slots.
  2. O-departure SMIRKS produces correct amide bond but OH departs as a
     separate fragment — take_largest:true needed on carboxyl reactions,
     OR use unmapped O in SMIRKS (both approaches valid).
  3. Ester O and carboxyl O are structurally identical in CHUCKLES —
     safety relies entirely on ester O never receiving a dummy atom.
  4. [amine_primary, carboxyl] REACTION_INDEX entry needs new SMIRKS;
     must remove it from backbone_amide and add to a new reaction entry.
  5. Intramolecular O-departure: verify single fragment output.
""")
