#!/usr/bin/env python3
"""
Validation of the pyPept reaction SMIRKS library.
Each reaction tested with [n*] dummy atoms as in CHUCKLES monomers.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

PASS = 0
FAIL = 0


def largest_fragment(mol):
    """Return the heaviest fragment of a disconnected molecule."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def run_smirks(smirks, *mols):
    rxn = AllChem.ReactionFromSmarts(smirks)
    if rxn is None:
        return None, f"invalid SMIRKS"
    n = rxn.GetNumReactantTemplates()
    if len(mols) != n:
        return None, f"{n} templates but {len(mols)} mols"
    products = rxn.RunReactants(tuple(mols))
    if not products:
        return None, "no products"
    results = []
    for ps in products:
        smi_set = []
        for mol in ps:
            try:
                Chem.SanitizeMol(mol)
                smi_set.append(Chem.MolToSmiles(mol))
            except Exception as e:
                smi_set.append(f"ERR({e})")
        results.append(smi_set)
    return results, None


def test(name, smirks_steps, reactant_smiles, expected_smiles,
         allow_multi=False, take_largest=False):
    global PASS, FAIL
    mols = [Chem.MolFromSmiles(s) for s in reactant_smiles]
    if any(m is None for m in mols):
        print(f"  SKIP  {name}: bad SMILES")
        return

    current_mol = None
    for step_i, smirks in enumerate(smirks_steps):
        if step_i == 0:
            results, err = run_smirks(smirks, *mols)
        else:
            results, err = run_smirks(smirks, current_mol)
        if err:
            print(f"  FAIL  {name} step {step_i+1}: {err}")
            FAIL += 1
            return
        unique = sorted({s for ps in results for s in ps if not s.startswith("ERR")})
        if not unique:
            print(f"  FAIL  {name} step {step_i+1}: all products failed sanitization")
            FAIL += 1
            return
        if not allow_multi and len(unique) > 1:
            print(f"  WARN  {name}: {len(unique)} regioisomers: {unique}")
        current_mol = Chem.MolFromSmiles(unique[0])

    if take_largest:
        current_mol = largest_fragment(current_mol)

    got = Chem.MolToSmiles(current_mol)
    exp_mol = Chem.MolFromSmiles(expected_smiles)
    exp = Chem.MolToSmiles(exp_mol) if exp_mol else expected_smiles

    if got == exp:
        print(f"  PASS  {name}")
        PASS += 1
    else:
        print(f"  FAIL  {name}")
        print(f"        expected: {exp}")
        print(f"        got:      {got}")
        FAIL += 1


# ── BACKBONE ────────────────────────────────────────────────────────
print("BACKBONE BONDS")

test("backbone_amide",
     ["[1*][N:2].[2*][C:4](=[O:5]) >> [N:2][C:4](=[O:5])"],
     ["[1*]NC", "CC(=O)[2*]"],
     "CNC(C)=O")

test("backbone_ester",
     ["[1*][O:2].[2*][C:4](=[O:5]) >> [O:2][C:4](=[O:5])"],
     ["[1*]OC", "CC(=O)[2*]"],
     "COC(C)=O")

# ── SIDECHAIN ───────────────────────────────────────────────────────
print("\nSIDECHAIN / CROSSLINK BONDS")

test("disulfide",
     ["[4*][S:2].[4*][S:4] >> [S:2][S:4]"],
     ["[4*]SC", "[4*]SC"],
     "CSSC")

test("thioester",
     ["[4*][S:2].[2*][C:4](=[O:5]) >> [S:2][C:4](=[O:5])"],
     ["[4*]SC", "CC(=O)[2*]"],
     "CSC(C)=O")

# thioether from alkyl halide — Cl departs, S-C formed
# [4*]SC = CH3-S-dummy; [4*]CCl = dummy-CH2-Cl
# Product: CH3-S + CH2 (Cl gone, dummy gone) → CH3-S-CH3 = CSC
test("thioether_halide",
     ["[4*][S:2].[4*][C:4][Cl,Br,I] >> [S:2][C:4]"],
     ["[4*]SC", "[4*]CCl"],
     "CSC")

# thiol-maleimide Michael addition — ring retained, C=C becomes C-C
test("thiol_maleimide",
     ["[4*][C:2]1=[C:3][C:4](=[O:5])[N:6][C:7](=[O:8])1.[4*][S:10] >> [C:2]1([S:10])[C:3][C:4](=[O:5])[N:6][C:7](=[O:8])1"],
     ["O=C1C([4*])=CC(=O)N1", "[4*]SC"],
     "O=C1CC(SC)C(=O)N1")

# NHS ester + amine — NHS ring departs as multi-atom leaving group
test("nhs_ester_amide",
     ["[4*][C:2](=[O:3])ON1C(=O)CCC1=O.[4*][N:5] >> [N:5][C:2](=[O:3])"],
     ["[4*]C(=O)ON1C(=O)CCC1=O", "[4*]NC"],
     "CNC=O")

# oxime ligation: aminooxy + aldehyde → C=N-O + (implicit H2O)
# [4*]ON = dummy-O-NH2; [4*]C=O = dummy-CHO
# O:2 bonded to N:3; C:5 loses O (aldehyde O unmapped → removed) and dummy
test("oxime_ligation",
     ["[4*][O:2][N:3].[4*][CH:5]=[O] >> [O:2][N:3]=[CH:5]"],
     ["[4*]ON", "[4*]C=O"],
     "[CH]=NO")

# hydrazone: hydrazide + aldehyde → C=N-N + (implicit H2O)
# [4*]NN = dummy-NH-NH2; [4*]C=O = dummy-CHO
test("hydrazone",
     ["[4*][N:2][N:3].[4*][C:4]=[O] >> [N:2][N:3]=[C:4]"],
     ["[4*]NN", "[4*]C=O"],
     "C=NN")

# ── CLICK CHEMISTRY ─────────────────────────────────────────────────
print("\nCLICK CHEMISTRY")

# CuAAC: terminal alkyne + azide → 1,4-substituted-1H-1,2,3-triazole
# Azide N atoms unmapped in reactant → unmapped in product (Gillingham approach)
# [4*] on inner alkyne C; azide alpha-C mapped; azide N atoms unmapped
test("CuAAC_1_4_triazole",
     ["[4*][C:2]#[CH:3].[4*][C:5][N]=[N+]=[N-] >> [C:5][n]1[n][n][c:3][c:2]1"],
     ["[4*]C#C", "[4*]CN=[N+]=[N-]"],
     "Cn1ccnn1",
     allow_multi=True)

# IEDDA two-step (simplified diazadiene model — real tetrazine+TCO needs scaffold-specific SMARTS)
# Step 1: [4+2] → bicyclic intermediate
# Step 2: retro-[4+2] N2 expulsion → diene + N#N
# take_largest=True: filter out N#N byproduct
test("IEDDA_2step_simplified",
     ["[C:1]=[N:2]-[N:3]=[C:4].[C:5]=[C:6] >> [C:1]1[N:2][N:3][C:4][C:6][C:5]1",
      "[C:1]1[N:2][N:3][C:4][C:6][C:5]1 >> [C:1]=[C:5][C:6]=[C:4].[N:2]#[N:3]"],
     ["C=NN=C", "C=C"],
     "C=CC=C",
     allow_multi=True, take_largest=True)

# ── RESTORATION (unbound sites) ─────────────────────────────────────
print("\nUNBOUND SITE RESTORATION")

test("restore_R1_H",
     ["[1*][N:2] >> [H][N:2]"],
     ["[1*][NH]C"],
     "CN")

test("restore_R2_OH",
     ["[2*][C:2](=[O:3]) >> O[C:2](=[O:3])"],
     ["CC(=O)[2*]"],
     "CC(=O)O")

test("restore_R3_H",
     ["[3*][N:2] >> [H][N:2]"],
     ["[3*][NH]CC"],
     "NCC")

test("restore_R4_thiol",
     ["[4*][S:2] >> [H][S:2]"],
     ["[4*]SC"],
     "CS")

test("restore_R4_amine",
     ["[4*][N:2] >> [H][N:2]"],
     ["[4*][NH]C"],
     "CN")

# ── SUMMARY ─────────────────────────────────────────────────────────
print(f"\nPASS {PASS}  FAIL {FAIL}  TOTAL {PASS+FAIL}")
sys.exit(0 if FAIL == 0 else 1)
