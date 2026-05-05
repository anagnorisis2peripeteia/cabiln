#!/usr/bin/env python3
"""Batch 2 monomer additions: Phe/Trp/Tyr variants, α-methyl AAs, β³-AAs,
Pro analogs, linkers, fatty acid caps, Glu/Asp variants.

Run from repo root:
    python tools/add_monomers_batch2.py

Monomers already in the SDF are silently skipped.
pre_activate is used for standard α-AAs; manual CHUCKLES are supplied for
β³-amino acids and cap-type monomers that pre_activate cannot parse.
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# ── Entry format ────────────────────────────────────────────────────────────
# (smiles_or_chuckles, abbr, name, m_type, m_subtype,
#  use_pre_activate,   # True = run pre_activate; False = smiles IS the CHUCKLES
#  m_Rgroups_override, m_chem_types_override)   # both None unless use_pre_activate=False
#
# For cap/linker monomers written as CHUCKLES directly:
#   m_Rgroups_override: "R1,R2,R3,R4,R5,R6" (None = not present)
#   m_chem_types_override: "slot:type,..." string

NEW_MONOMERS = [

    # ── Phe ring variants ────────────────────────────────────────────────────
    ("N[C@@H](Cc1ccc(C(F)(F)F)cc1)C(=O)O",  "Phe_4CF3",  "L-4-(Trifluoromethyl)phenylalanine", "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(C#N)cc1)C(=O)O",       "Phe_4CN",   "L-4-Cyanophenylalanine",             "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(OC)cc1)C(=O)O",        "Phe_4OMe",  "L-4-Methoxyphenylalanine",           "aa","modified", True, None, None),
    ("N[C@@H](Cc1cccc(OC)c1)C(=O)O",        "Phe_3OMe",  "L-3-Methoxyphenylalanine",           "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccccc1OC)C(=O)O",          "Phe_2OMe",  "L-2-Methoxyphenylalanine",           "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccccc1Cl)C(=O)O",          "Phe_2Cl",   "L-2-Chlorophenylalanine",            "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccccc1Br)C(=O)O",          "Phe_2Br",   "L-2-Bromophenylalanine",             "aa","modified", True, None, None),
    ("N[C@@H](Cc1cccc(Br)c1)C(=O)O",        "Phe_3Br",   "L-3-Bromophenylalanine",             "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(F)c(F)c1)C(=O)O",      "Phe_34diF", "L-3,4-Difluorophenylalanine",        "aa","modified", True, None, None),
    ("N[C@@H](Cc1cc(F)cc(F)c1)C(=O)O",      "Phe_35diF", "L-3,5-Difluorophenylalanine",        "aa","modified", True, None, None),
    ("N[C@@H](Cc1c(F)cccc1F)C(=O)O",        "Phe_26diF", "L-2,6-Difluorophenylalanine",        "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(F)cc1F)C(=O)O",        "Phe_24diF", "L-2,4-Difluorophenylalanine",        "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(C(C)(C)C)cc1)C(=O)O",  "Phe_4tBu",  "L-4-tert-Butylphenylalanine",        "aa","modified", True, None, None),
    ("N[C@@H](Cc1c(F)c(F)c(F)c(F)c1F)C(=O)O", "Phe_F5", "L-Pentafluorophenylalanine",          "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(C(C)=O)cc1)C(=O)O",    "Phe_4Ac",   "L-4-Acetylphenylalanine",            "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccccc1[N+](=O)[O-])C(=O)O","Phe_2NO2",  "L-2-Nitrophenylalanine",             "aa","modified", True, None, None),
    ("N[C@@H](Cc1cccc([N+](=O)[O-])c1)C(=O)O","Phe_3NO2","L-3-Nitrophenylalanine",             "aa","modified", True, None, None),
    ("N[C@@H](Cc1cccc(O)c1)C(=O)O",         "Phe_3OH",   "L-3-Hydroxyphenylalanine",           "aa","modified", True, None, None),
    ("N[C@@H](CCc1ccccc1)C(=O)O",           "hPhe2",     "L-Homophenylalanine (2-amino-4-phenylbutanoic acid)", "aa","modified", True, None, None),

    # ── Trp ring variants ───────────────────────────────────────────────────
    ("N[C@@H](Cc1c[nH]c2cc(F)ccc12)C(=O)O", "Trp_5F",   "L-5-Fluorotryptophan",               "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2cc(C)ccc12)C(=O)O", "Trp_5Me",  "L-5-Methyltryptophan",               "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2cc(Br)ccc12)C(=O)O","Trp_5Br",  "L-5-Bromotryptophan",                "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2cc(OC)ccc12)C(=O)O","Trp_5OMe", "L-5-Methoxytryptophan",              "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O", "Trp_6F",   "L-6-Fluorotryptophan",               "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2ncccc12)C(=O)O",    "Trp_7Aza", "L-7-Azatryptophan",                  "aa","modified", True, None, None),
    ("N[C@@H](Cc1cn(C)c2ccccc12)C(=O)O",    "Trp_1Me",  "L-1-Methyltryptophan",               "aa","modified", True, None, None),
    ("N[C@@H](Cc1c[nH]c2c(F)cccc12)C(=O)O", "Trp_4F",   "L-4-Fluorotryptophan",               "aa","modified", True, None, None),
    ("N[C@@H](CCc1c[nH]c2ccccc12)C(=O)O",   "hTrp",     "L-Homotryptophan",                   "aa","modified", True, None, None),

    # ── Tyr ring variants ───────────────────────────────────────────────────
    ("N[C@@H](Cc1ccc(O)c(F)c1)C(=O)O",      "Tyr_3F",   "L-3-Fluorotyrosine",                 "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(O)c(Cl)c1)C(=O)O",     "Tyr_3Cl",  "L-3-Chlorotyrosine",                 "aa","modified", True, None, None),
    ("N[C@@H](Cc1ccc(O)c(Br)c1)C(=O)O",     "Tyr_3Br",  "L-3-Bromotyrosine",                  "aa","modified", True, None, None),
    ("N[C@@H](Cc1c(C)cc(O)cc1C)C(=O)O",     "Tyr_35diMe","L-3,5-Dimethyltyrosine",             "aa","modified", True, None, None),

    # ── Alpha-methyl amino acids ─────────────────────────────────────────────
    ("NC(C)(C(C)C)C(=O)O",                   "aMeVal",   "alpha-Methylvaline (alpha-Me-Val)",   "aa","modified", True, None, None),
    ("NC(C)(Cc1c[nH]c2ccccc12)C(=O)O",      "aMeTrp",   "alpha-Methyltryptophan",              "aa","modified", True, None, None),
    ("NC(C)(CCC(=O)O)C(=O)O",               "aMeGlu",   "alpha-Methylglutamic acid",           "aa","modified", True, None, None),
    ("NC(C)(CC(=O)O)C(=O)O",                "aMeAsp",   "alpha-Methylaspartic acid",           "aa","modified", True, None, None),
    ("NC(C)(CCCCN)C(=O)O",                  "aMeLys",   "alpha-Methyllysine",                  "aa","modified", True, None, None),
    ("NC(C)(CCC(=O)N)C(=O)O",               "aMeGln",   "alpha-Methylglutamine",               "aa","modified", True, None, None),
    ("NC(C)(CC(=O)N)C(=O)O",                "aMeAsn",   "alpha-Methylasparagine",              "aa","modified", True, None, None),
    ("NC(C)(CO)C(=O)O",                     "aMeSer",   "alpha-Methylserine",                  "aa","modified", True, None, None),
    ("NC(C)(CS)C(=O)O",                     "aMeCys",   "alpha-Methylcysteine",                "aa","modified", True, None, None),
    ("NC(C)(Cc1cnc[nH]1)C(=O)O",            "aMeHis",   "alpha-Methylhistidine",               "aa","modified", True, None, None),
    ("NC(C)([C@@H](C)CC)C(=O)O",            "aMeIle",   "alpha-Methylisoleucine",              "aa","modified", True, None, None),
    ("NC(C)(CCCNC(=N)N)C(=O)O",             "aMeArg",   "alpha-Methylarginine",                "aa","modified", True, None, None),

    # ── β³-Amino acids (manual CHUCKLES — pre_activate doesn't find β-backbone) ─
    # CHUCKLES pattern: [1*]N([3*])[C@@H](sidechain)CC([2*])=O
    # m_Rgroups: [H],[OH],[H],<sidechain_R or None>,None,None
    # Note: b3hLys and b3hGlu have sidechain R4/R3 — handled via their CHUCKLES

    ("[1*]N([3*])[C@@H](Cc1ccccc1)CC([2*])=O",          "b3hPhe", "beta3-Homophenylalanine",      "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](CC(C)C)CC([2*])=O",             "b3hLeu", "beta3-Homoleucine",            "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](C(C)C)CC([2*])=O",              "b3hVal", "beta3-Homovaline",             "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](Cc1c[nH]c2ccccc12)CC([2*])=O",  "b3hTrp", "beta3-Homotryptophan",        "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](Cc1ccc(O)cc1)CC([2*])=O",       "b3hTyr", "beta3-Homotyrosine",          "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](CO)CC([2*])=O",                 "b3hSer", "beta3-Homoserine",             "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](CCCCN[4*])CC([2*])=O",          "b3hLys", "beta3-Homolysine",            "aa","modified", False,
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary"),

    ("[1*]N([3*])[C@@H](CCC(=O)O[4*])CC([2*])=O",       "b3hGlu", "beta3-Homoglutamic acid",     "aa","modified", False,
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@@H](CC(=O)O[4*])CC([2*])=O",        "b3hAsp", "beta3-Homoaspartic acid",     "aa","modified", False,
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@@H](C)CC([2*])=O",                  "b3hAla", "beta3-Homoalanine",           "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@@H](CCCC)CC([2*])=O",               "b3hNle", "beta3-Homonorleucine",        "aa","modified", False,
     "[H],[OH],[H],None,None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    # ── Pro / cyclic analogs ─────────────────────────────────────────────────
    ("OC(=O)[C@@H]1C[C@H](F)CN1",           "Pro_4F_t", "trans-4-Fluoroproline (2S,4R)",       "aa","modified", True, None, None),
    ("OC(=O)[C@@H]1C[C@@H](F)CN1",          "Pro_4F_c", "cis-4-Fluoroproline (2S,4S)",         "aa","modified", True, None, None),
    ("OC(=O)[C@@H]1Cc2ccccc2N1",            "Ind",      "L-Indoline-2-carboxylic acid",        "aa","modified", True, None, None),

    # ── Glu / Asp side-chain protected variants ──────────────────────────────
    ("N[C@@H](CCC(=O)OC)C(=O)O",            "Glu_OMe",  "L-Glu gamma-methyl ester",            "aa","modified", True, None, None),
    ("N[C@@H](CCC(=O)OC(C)(C)C)C(=O)O",     "Glu_OtBu", "L-Glu gamma-tert-butyl ester",       "aa","modified", True, None, None),
    ("N[C@@H](CCC(=O)OCc1ccccc1)C(=O)O",    "Glu_OBn",  "L-Glu gamma-benzyl ester",           "aa","modified", True, None, None),
    ("N[C@@H](CC(=O)OC(C)(C)C)C(=O)O",      "Asp_OtBu", "L-Asp beta-tert-butyl ester",        "aa","modified", True, None, None),
    ("N[C@@H](CC(=O)OCc1ccccc1)C(=O)O",     "Asp_OBn",  "L-Asp beta-benzyl ester",            "aa","modified", True, None, None),

    # ── Linker / spacer monomers ─────────────────────────────────────────────
    ("NCCCC(=O)O",                           "GABA",    "gamma-Aminobutyric acid (4-aminobutyric acid)", "linker","linker", True, None, None),
    ("NCCCCC(=O)O",                          "5Ava",    "5-Aminopentanoic acid (delta-aminovaleric acid)", "linker","linker", True, None, None),
    ("NCCCCCC(=O)O",                         "6Ahx",    "6-Aminohexanoic acid (epsilon-aminocaproic acid)", "linker","linker", True, None, None),
    ("NCCCCCCC(=O)O",                        "7Ahp",    "7-Aminoheptanoic acid",               "linker","linker", True, None, None),
    ("NCCCCCCCC(=O)O",                       "8Ado",    "8-Aminooctanoic acid",                "linker","linker", True, None, None),
    ("NCCCCCCCCCCC(=O)O",                    "11Aun",   "11-Aminoundecanoic acid",             "linker","linker", True, None, None),
    ("NCCOCC(=O)O",                          "OEG1",    "5-Amino-3-oxapentanoic acid (mini-PEG 1-unit)", "linker","linker", True, None, None),
    ("NCCOCCOCCOCC(=O)O",                    "OEG3",    "11-Amino-3,6,9-trioxaundecanoic acid (mini-PEG 3-unit)", "linker","linker", True, None, None),
    ("NCCOCCOCCOCCOCC(=O)O",                 "OEG4",    "14-Amino-3,6,9,12-tetraoxatetradecanoic acid (mini-PEG 4-unit)", "linker","linker", True, None, None),

    # ── Fatty acid / diacid caps (manual CHUCKLES — no amine for pre_activate) ─
    # These are diacid caps: R2=carboxyl (bonds to linker amine), free COOH at other end
    ("[2*]OC(=O)CCCCCCCCCCCC(=O)O",         "C14dAcid", "Tetradecanedioic acid cap",           "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCC(=O)O",       "C16dAcid", "Hexadecanedioic acid cap",            "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCCCCCCCC(=O)O", "C22dAcid", "Docosanedioic acid cap",              "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    # Monoacid fatty acid caps: palmitic, myristic, stearic, oleic
    ("[2*]OC(=O)CCCCCCCCCCCCCC",             "Myr",     "Myristoyl cap (C14 monoacid)",        "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCCCC",           "Pal",     "Palmitoyl cap (C16 monoacid)",        "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCCCCCC",         "Ste",     "Stearoyl cap (C18 monoacid)",         "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCC/C=C\\CCCCCCCC",     "Ole",     "Oleoyl cap (C18:1 monoacid)",         "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCC(=O)O",       "C14FA",   "Tetradecanoic diacid cap (C14)",      "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCCCC(=O)O",     "C16FA",   "Hexadecanoic diacid cap (C16)",       "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),

    ("[2*]OC(=O)CCCCCCCCCCCCCCCCCCCCCC(=O)O","C22FA",  "Docosanoic diacid cap (C22)",         "cap","cap", False,
     "None,[OH],None,None,None,None", "2:backbone_c"),
]


def _format_chem_types(chem_types_dict: dict) -> str:
    return ",".join(f"{k}:{v}" for k, v in sorted(chem_types_dict.items()))


def _format_rgroups(leaving_dict: dict) -> str:
    slots = ["None"] * 6
    for slot, lg in leaving_dict.items():
        idx = int(slot) - 1
        if 0 <= idx < 6:
            slots[idx] = lg
    return ",".join(slots)


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor

    rdDepictor.SetPreferCoordGen(True)

    # Load existing abbreviations
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    existing = set()
    for m in suppl:
        if m:
            existing.add(m.GetPropsAsDict().get("m_abbr", ""))
    print(f"Existing monomers: {len(existing)}")

    added = 0
    skipped = 0
    failed = []

    with open(SDF_PATH, "a", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for entry in NEW_MONOMERS:
            smiles_or_chuckles, abbr, name, mtype, msubtype, use_pre_activate, rgroups_override, chem_types_override = entry

            if abbr in existing:
                skipped += 1
                continue

            try:
                if use_pre_activate:
                    from pyPept.interfaces.monomer_pipeline import pre_activate
                    result = pre_activate(smiles_or_chuckles)
                    mol = Chem.MolFromSmiles(result.chuckles)
                    rgroups_str = _format_rgroups(result.leaving)
                    chem_types_str = _format_chem_types(result.chem_types)
                else:
                    # smiles_or_chuckles IS the CHUCKLES; use directly
                    mol = Chem.MolFromSmiles(smiles_or_chuckles)
                    rgroups_str = rgroups_override
                    chem_types_str = chem_types_override

                if mol is None:
                    failed.append((abbr, "MolFromSmiles returned None"))
                    continue

                rdDepictor.Compute2DCoords(mol)
                mol.SetProp("m_abbr",       abbr)
                mol.SetProp("symbol",       abbr)
                mol.SetProp("m_name",       name)
                mol.SetProp("m_type",       mtype)
                mol.SetProp("m_subtype",    msubtype)
                mol.SetProp("m_Rgroups",    rgroups_str)
                mol.SetProp("m_chem_types", chem_types_str)
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")

            except Exception as exc:
                failed.append((abbr, str(exc).split("\n")[0]))

        writer.close()

    print(f"\nDone. Added {added} monomers (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")


if __name__ == "__main__":
    main()
