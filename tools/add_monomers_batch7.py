#!/usr/bin/env python3
"""Batch 7: final push to reach 1000 monomers.

Adds ~35 entries (expecting ~31 new) to bring 969 -> 1000.

Run from repo root:
    python tools/add_monomers_batch7.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

_AA_RG = "[H],[OH],[H],None,None,None"
_AA_CT = "1:backbone_n,2:backbone_c,3:backbone_n_mod"
_NME_RG = "[H],[OH],None,None,None,None"
_NME_CT = "1:backbone_n,2:backbone_c"
_NME1_RG = "[H],[OH],[H],None,None,None"
_NME1_CT = "1:backbone_n,2:backbone_c"
_AME_RG = "[H],[OH],[H],None,None,None"
_AME_CT = "1:backbone_n,2:backbone_c,3:backbone_n_mod"
_PRO_RG = "[H],[OH],None,None,None,None"
_PRO_CT = "1:backbone_n,2:backbone_c"
_CAP_RG = "None,[OH],None,None,None,None"
_CAP_CT = "2:backbone_c"

BATCH7 = [
    # ─── D-bAbu (beta-amino butyric acid, D-form of bAbu added in batch6) ───
    ("[1*]N([3*])[C@H](C)CC([2*])=O",
     "D_bAbu", "D-3-aminobutyric acid (D-beta-aminobutyric acid)", "aa", "modified", _AA_RG, _AA_CT),

    # ─── homo-Asp and homo-Glu (L and D) ──────────────────────────────────────
    ("[1*]N([3*])[C@@H](CCC(=O)O[4*])C([2*])=O",
     "hAsp", "L-homoaspartate (2-amino-pentanedioic 1-acid)", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@H](CCC(=O)O[4*])C([2*])=O",
     "D_hAsp", "D-homoaspartate", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@@H](CCCC(=O)O[4*])C([2*])=O",
     "hGlu", "L-homoglutamate (2-aminoadipic acid)", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@H](CCCC(=O)O[4*])C([2*])=O",
     "D_hGlu", "D-homoglutamate", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    # ─── homo-Pro = pipecolic acid variants (L/D already exist as Pip/D_Pip) ─
    # hTrp-5F (homo-5-fluorotryptophan, extra CH2 vs Trp_5F)
    ("[1*]N([3*])[C@@H](CCc1cn([4*])c2cc(F)ccc12)C([2*])=O",
     "hTrp_5F", "L-homo-5-fluorotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@H](CCc1cn([4*])c2cc(F)ccc12)C([2*])=O",
     "D_hTrp_5F", "D-homo-5-fluorotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    # ─── New L/D aMePhe 4-NO2 and 4-NH2 ─────────────────────────────────────
    ("[1*]N([3*])[C@@](C)(Cc1ccc([N+](=O)[O-])cc1)C([2*])=O",
     "aMePhe_4NO2", "L-alpha-methyl-4-nitrophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc([N+](=O)[O-])cc1)C([2*])=O",
     "D_aMePhe_4NO2", "D-alpha-methyl-4-nitrophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(N)cc1)C([2*])=O",
     "aMePhe_4NH2", "L-alpha-methyl-4-aminophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(N)cc1)C([2*])=O",
     "D_aMePhe_4NH2", "D-alpha-methyl-4-aminophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    # ─── b3h-Phe ortho/meta fluorinated ──────────────────────────────────────
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccccc1F",
     "b3hPhe_2F", "L-beta3-homo-2-fluorophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccccc1F",
     "D_b3hPhe_2F", "D-beta3-homo-2-fluorophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1cccc(F)c1",
     "b3hPhe_3F", "L-beta3-homo-3-fluorophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cccc(F)c1",
     "D_b3hPhe_3F", "D-beta3-homo-3-fluorophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # ─── NMe-Phe 3,5-diF and 3,4-diCl ──────────────────────────────────────
    ("[1*]N(C)[C@@H](Cc1cc(F)cc(F)c1)C([2*])=O",
     "mePhe_35diF", "NMe-L-3,5-difluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1cc(F)cc(F)c1)C([2*])=O",
     "D_mePhe_35diF", "NMe-D-3,5-difluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(Cl)c(Cl)c1)C([2*])=O",
     "mePhe_34diCl", "NMe-L-3,4-dichlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(Cl)c(Cl)c1)C([2*])=O",
     "D_mePhe_34diCl", "NMe-D-3,4-dichlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # ─── More aMeTrp and related ─────────────────────────────────────────────
    ("[1*]N([3*])[C@@](C)(Cc1cn([4*])c2cc(Br)ccc12)C([2*])=O",
     "aMeTrp_5Br", "L-alpha-methyl-5-bromotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@](C)(Cc1cn([4*])c2cc(Br)ccc12)C([2*])=O",
     "D_aMeTrp_5Br", "D-alpha-methyl-5-bromotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    # ─── Additional useful AAs ───────────────────────────────────────────────
    # Ala_3F (beta-fluoroalanine)
    ("[1*]N([3*])[C@@H](CF)C([2*])=O",
     "Ala_3F", "L-beta-fluoroalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CF)C([2*])=O",
     "D_Ala_3F", "D-beta-fluoroalanine", "aa", "modified", _AA_RG, _AA_CT),

    # Kynurenine (Kyn) - aromatic AA from Trp metabolism
    ("[1*]N([3*])[C@@H](CC(=O)c1ccccc1N)C([2*])=O",
     "Kyn", "L-kynurenine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC(=O)c1ccccc1N)C([2*])=O",
     "D_Kyn", "D-kynurenine", "aa", "modified", _AA_RG, _AA_CT),

    # Hydroxyproline 4-OBn (protected)
    ("[1*]N1C[C@@H](OCc2ccccc2)C[C@H]1C([2*])=O",
     "Pro_4OBn", "L-trans-4-benzyloxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO_RG, _PRO_CT),

    ("[1*]N1C[C@@H](OCc2ccccc2)C[C@@H]1C([2*])=O",
     "D_Pro_4OBn", "D-trans-4-benzyloxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO_RG, _PRO_CT),

    # Thienyl alanine (2-thienyl and 3-thienyl)
    ("[1*]N([3*])[C@@H](Cc1cccs1)C([2*])=O",
     "Thi2Ala", "L-2-thienylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1cccs1)C([2*])=O",
     "D_Thi2Ala", "D-2-thienylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](Cc1ccsc1)C([2*])=O",
     "Thi3Ala", "L-3-thienylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1ccsc1)C([2*])=O",
     "D_Thi3Ala", "D-3-thienylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # Furylalanine (2-furyl)
    ("[1*]N([3*])[C@@H](Cc1ccco1)C([2*])=O",
     "Fur2Ala", "L-2-furylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1ccco1)C([2*])=O",
     "D_Fur2Ala", "D-2-furylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # Aminocaproic acid (Ahx) as linker
    ("[1*]N([3*])CCCCC([2*])=O",
     "Ahx", "6-aminohexanoic acid (Ahx) linker", "linker", "linker",
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    # Trioxadecanoic acid (mini-PEG linker, OEG3 variant)
    ("[1*]N([3*])CCOCCOCCC([2*])=O",
     "OEG3", "8-amino-3,6-dioxaoctanoic acid (OEG3 PEG linker)", "linker", "linker",
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),
]


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor

    rdDepictor.SetPreferCoordGen(True)

    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    existing = set()
    for m in suppl:
        if m:
            existing.add(m.GetPropsAsDict().get("m_abbr", ""))
    print(f"Existing monomers: {len(existing)}")

    added = 0
    skipped = 0
    failed = []

    with open(str(SDF_PATH), "a", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for chuckles, abbr, name, mtype, msubtype, rg, ct in BATCH7:
            if abbr in existing:
                skipped += 1
                continue
            try:
                mol = Chem.MolFromSmiles(chuckles)
                if mol is None:
                    raise ValueError("MolFromSmiles returned None")
                rdDepictor.Compute2DCoords(mol)
                mol.SetProp("m_abbr", abbr)
                mol.SetProp("symbol", abbr)
                mol.SetProp("m_name", name)
                mol.SetProp("m_type", mtype)
                mol.SetProp("m_subtype", msubtype)
                mol.SetProp("m_Rgroups", rg)
                mol.SetProp("m_chem_types", ct)
                writer.write(mol)
                added += 1
            except Exception as exc:
                failed.append((abbr, str(exc)[:120]))

        writer.close()

    print(f"Done. Added {added} (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")

    s = list(Chem.SDMolSupplier(str(SDF_PATH), removeHs=False))
    valid = [m for m in s if m]
    print(f"\nFinal valid monomers: {len(valid)}")


if __name__ == "__main__":
    main()
