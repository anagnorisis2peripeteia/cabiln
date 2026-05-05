#!/usr/bin/env python3
"""Test reagent-form cap SMILES through the pipeline."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.interfaces.monomer_pipeline import pre_activate, find_cap_slots, infer_leaving_group
from rdkit import Chem

test_cases = [
    # (reagent SMILES, name, expected LG, expected slot)
    # --- Acyl chlorides (reagent form of acylation caps) ---
    ("CC(=O)Cl",                    "acetyl chloride (ac)",       "[Cl]", 2),
    ("CC(=O)O",                     "acetic acid (ac free form)", "[OH]", 2),
    ("ClC(=O)OCc1ccccc1",           "Cbz-Cl",                    "[Cl]", 2),
    ("CC(C)(C)OC(=O)Cl",            "Boc-Cl",                    "[Cl]", 2),
    ("O=C(Cl)Cc1ccccc1",            "phenylacetyl chloride (BnCap)", "[Cl]", 2),
    ("O=C(Cl)CBr",                  "bromoacetyl chloride (BrAc)","[Cl]", 2),
    ("O=C(Cl)CCl",                  "chloroacetyl chloride (ClAc)","[Cl]", 2),
    ("O=C(Cl)C(F)(F)F",             "trifluoroacetyl Cl (CF3Ac)", "[Cl]", 2),
    ("O=C(Cl)c1ccccc1",             "benzoyl chloride (Bz)",      "[Cl]", 2),
    ("O=C(Cl)C(C)(C)C",             "pivaloyl chloride (PivCap)", "[Cl]", 2),

    # --- Sulfonyl chlorides ---
    ("Cc1ccc(S(=O)(=O)Cl)cc1",      "tosyl chloride (Tos)",       "[Cl]", 2),
    ("CS(=O)(=O)Cl",                 "mesyl chloride (MsO)",       "[Cl]", 2),

    # --- Alkyl halides ---
    ("CI",                           "methyl iodide (Me_)",        "[I]",  2),
    ("CCBr",                         "ethyl bromide (Et_)",        "[Br]", 2),
    ("BrCc1ccccc1",                  "benzyl bromide (Bn_)",       "[Br]", 2),
    ("ClC(c1ccccc1)(c1ccccc1)c1ccccc1", "trityl chloride (trt)",  "[Cl]", 2),

    # --- Nucleophilic C-terminal caps ---
    ("N",                            "ammonia (am cap)",           "[H]",  1),
    ("CN",                           "methylamine (NHMe cap)",     "[H]",  1),
    ("CCN",                          "ethylamine (NHEt cap)",      "[H]",  1),
    ("NCc1ccccc1",                   "benzylamine (_NHBn)",        "[H]",  1),
    ("CO",                           "methanol (_OMe cap)",        "[H]",  1),
    ("CCO",                          "ethanol (OEt cap)",          "[H]",  1),
    ("OCc1ccccc1",                   "benzyl alcohol (_OBn)",      "[H]",  1),

    # --- Free-form acylation caps (HATU coupling) ---
    ("OC(=O)CCCCCCCCCCCCC",          "myristic acid (Myr)",        "[OH]", 2),
    ("OC(=O)CCCCC(=O)O",             "adipic acid (Adp)",          "[OH]", 2),
    ("OC(=O)CCC(=O)O",               "succinic acid (Succ)",       "[OH]", 2),
]

print("=" * 80)
print("REAGENT-FORM CAP DETECTION TEST")
print("=" * 80)

passed = 0
failed = 0

for smi, name, expected_lg, expected_slot in test_cases:
    print(f"\n--- {name}: {smi} ---")

    try:
        result = pre_activate(smi)
        if expected_slot in result.leaving and result.leaving[expected_slot] == expected_lg:
            print(f"  PASS: chuckles={result.chuckles}  LG={result.leaving}")
            passed += 1
        else:
            print(f"  FAIL: expected slot {expected_slot} LG={expected_lg}, "
                  f"got leaving={result.leaving}")
            print(f"        chuckles={result.chuckles}")
            failed += 1
    except Exception as e:
        print(f"  FAIL (exception): {e}")
        failed += 1

print(f"\n{'=' * 80}")
print(f"Results: {passed} passed, {failed} failed out of {len(test_cases)}")
print(f"{'=' * 80}")
