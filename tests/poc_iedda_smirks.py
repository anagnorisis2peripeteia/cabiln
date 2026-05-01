#!/usr/bin/env python3
"""
Proof-of-concept: IEDDA two-step SMIRKS via RDKit RunReactants.

Tests three things in order:
  1. Disconnected product (N2 expulsion) — the key feasibility question
  2. Ring-forming [4+2] Diels-Alder
  3. Chained two-step: [4+2] then retro-[4+2] N2 expulsion

Uses simplified models (no strained rings) to isolate RDKit capability
from SMARTS-writing complexity.
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def run(label, smirks, *reactant_smiles):
    rxn = AllChem.ReactionFromSmarts(smirks)
    assert rxn is not None, f"{label}: bad SMIRKS"
    mols = tuple(Chem.MolFromSmiles(s) for s in reactant_smiles)
    assert all(m is not None for m in mols), f"{label}: bad SMILES"
    products = rxn.RunReactants(mols)
    if not products:
        print(f"  FAIL {label}: no products")
        return []
    out = []
    for product_set in products:
        smi_set = []
        for mol in product_set:
            try:
                Chem.SanitizeMol(mol)
                smi_set.append(Chem.MolToSmiles(mol))
            except Exception as e:
                smi_set.append(f"SANITIZE_ERROR({e})")
        out.append(smi_set)
    print(f"  OK   {label}")
    for i, s in enumerate(out):
        print(f"       product set {i}: {s}")
    return out


print("=" * 60)
print("TEST 1: Disconnected product — N2 expulsion from azo compound")
print("=" * 60)
# Reactant: cis-azomethane CH3-N=N-CH3 (azo, thermal → 2 CH3• + N2, but
# we model the concerted retro as: break both C-N, form N≡N, carbons become
# radicals — use alkyl chain model where C-C forms instead)
# Simpler: cyclic azo → N2 + alkene
# 3,4-diaza-3-cyclopentene: C1CC=NN1  ->  butadiene + N2
smirks_n2 = '[C:1][N:2]=[N:3][C:4] >> [C:1]=[C:4].[N:2]#[N:3]'
run("N2 expulsion (acyclic azo, C=C + N2)",
    smirks_n2, 'C/N=N/C')

# Cyclic: break both C-N bonds in a ring to get diene + N2
smirks_n2_ring = '[C:1]1[C:2][N:3]=[N:4][C:5]1 >> [C:1]=[C:2].[N:3]#[N:4].[C:5]'
# Simpler version — just the N2 leaving from a ring without worrying about all atoms
smirks_n2_ring2 = '[C:1]-[N:2]=[N:3]-[C:4] >> [C:1]=[C:4].[N:2]#[N:3]'
run("N2 expulsion (acyclic azo, same atoms)", smirks_n2_ring2, 'CC/N=N/CC')


print()
print("=" * 60)
print("TEST 2: Ring-forming [4+2] Diels-Alder")
print("=" * 60)
# All-carbon DA: butadiene + ethylene -> cyclohexene
smirks_da = '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6] >> [C:1]1[C:2]=[C:3][C:4][C:6][C:5]1'
run("All-carbon DA (butadiene + ethylene -> cyclohexene)",
    smirks_da, 'C=CC=C', 'C=C')

# Hetero DA: C=N-N=C diene (tetrazine fragment) + ethylene
# This is the key pattern for IEDDA step 1
smirks_hda = '[C:1]=[N:2][N:3]=[C:4].[C:5]=[C:6] >> [C:1]1[N:2][N:3][C:4][C:6][C:5]1'
run("Hetero DA (C=N-N=C + ethylene) — IEDDA step 1 model",
    smirks_hda, 'C(/N=N/C)=C\\C', 'C=C')

# Even simpler: 1-aza-1,3-butadiene + ethylene
smirks_azada = '[C:1]=[N:2]-[C:3]=[C:4].[C:5]=[C:6] >> [C:1]1[N:2][C:3]=[C:4][C:6][C:5]1'
run("Aza-DA (C=N-C=C + ethylene)",
    smirks_azada, 'C=NC=C', 'C=C')


print()
print("=" * 60)
print("TEST 3: Two-step IEDDA — [4+2] then N2 expulsion")
print("=" * 60)
# Minimal model:
# Reactant 1: a diazadiene (mimics tetrazine): C=N-N=C
# Reactant 2: ethylene
# Step 1: [4+2] -> bicyclic intermediate with N-N bridge
# Step 2: retro [4+2] expels N2, forms new alkene

# Step 1 SMIRKS: C=N-N=C + C=C -> six-membered ring with N-N still present
smirks_step1 = '[C:1]=[N:2]-[N:3]=[C:4].[C:5]=[C:6] >> [C:1]1[N:2][N:3][C:4][C:6][C:5]1'

# Step 2 SMIRKS: break the N-N bridge out of the ring -> dihydropyridazine + N2
# In the 6-membered ring [C1-N2-N3-C4-C6-C5]:
# Break C1-N2 and C4-N3, form N2#N3, C1=C4... wait, C1 and C4 aren't adjacent
# Actually: break C1-N2 and C4-N3 bonds, N2-N3 becomes N2 gas
# Product ring: C1-C5-C6-C4 (4 atoms) — that's only 4 atoms left, forms butadiene fragment
# Let's check: if we remove N2-N3 from the 6-membered ring [C1-N2-N3-C4-C6-C5],
# we get fragments C1 and C4, connected through C5-C6
# So product would be: C1=C5-C6=C4 (a 1,3-diene) + N2
smirks_step2 = '[C:1]1[N:2][N:3][C:4][C:6][C:5]1 >> [C:1]=[C:5][C:6]=[C:4].[N:2]#[N:3]'

print("  Step 1 SMIRKS:", smirks_step1)
print("  Step 2 SMIRKS:", smirks_step2)

tetrazine_model = 'C(/N=N/C)=C'   # simplified: trans-C(=C)-N=N-C fragment
# Actually need a proper diazadiene: CH2=N-N=CH2 (diazomethylene model)
diazadiene = 'C=NN=C'
dienophile = 'C=C'

step1_products = run("IEDDA step 1", smirks_step1, diazadiene, dienophile)
if step1_products:
    intermediate_smi = step1_products[0][0]
    print(f"  Intermediate: {intermediate_smi}")
    step2_products = run("IEDDA step 2 (N2 expulsion)", smirks_step2, intermediate_smi)
    if step2_products:
        print(f"  Final products: {step2_products[0]}")
    else:
        print("  Step 2 FAILED — trying alternative SMIRKS")
        # Alternative: match N-N within ring more explicitly
        smirks_step2b = '[C:1]1-[N:2]-[N:3]-[C:4]-[C:6]-[C:5]-1 >> [C:1]=[C:5].[C:4]=[C:6].[N:2]#[N:3]'
        run("IEDDA step 2 alt", smirks_step2b, intermediate_smi)


print()
print("=" * 60)
print("TEST 4: RDKit handles N#N as product?")
print("=" * 60)
# Sanity check: can we even produce N#N in a SMIRKS product?
smirks_n2_simple = '[N:1]-[N:2] >> [N:1]#[N:2]'
run("N=N -> N#N bond order change", smirks_n2_simple, 'N=N')

smirks_n2_disconnect = '[C:1]-[N:2]-[N:3]-[C:4] >> [C:1].[N:2]#[N:3].[C:4]'
run("C-N-N-C -> C + N#N + C (disconnected)", smirks_n2_disconnect, 'CNNC')
