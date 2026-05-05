#!/usr/bin/env python3
"""Batch 4: D-forms of batch2 L-AAs, more Trp/Pro ring variants, unusual AAs, PEG series.

Run from repo root:
    python tools/add_monomers_batch4.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# (smiles, abbr, name, m_type, m_subtype)
# All use pre_activate — standard alpha-AAs with free NH2 and COOH.
NEW_MONOMERS = [

    # ── D-Phe ring variants ───────────────────────────────────────────────────
    ("N[C@H](Cc1ccc(C(F)(F)F)cc1)C(=O)O",    "D_Phe_4CF3",  "D-4-(Trifluoromethyl)phenylalanine",    "aa","modified"),
    ("N[C@H](Cc1ccc(C#N)cc1)C(=O)O",          "D_Phe_4CN",   "D-4-Cyanophenylalanine",                "aa","modified"),
    ("N[C@H](Cc1ccc(OC)cc1)C(=O)O",           "D_Phe_4OMe",  "D-4-Methoxyphenylalanine",              "aa","modified"),
    ("N[C@H](Cc1cccc(OC)c1)C(=O)O",           "D_Phe_3OMe",  "D-3-Methoxyphenylalanine",              "aa","modified"),
    ("N[C@H](Cc1ccccc1OC)C(=O)O",             "D_Phe_2OMe",  "D-2-Methoxyphenylalanine",              "aa","modified"),
    ("N[C@H](Cc1ccccc1Cl)C(=O)O",             "D_Phe_2Cl",   "D-2-Chlorophenylalanine",               "aa","modified"),
    ("N[C@H](Cc1ccccc1Br)C(=O)O",             "D_Phe_2Br",   "D-2-Bromophenylalanine",                "aa","modified"),
    ("N[C@H](Cc1cccc(Br)c1)C(=O)O",           "D_Phe_3Br",   "D-3-Bromophenylalanine",                "aa","modified"),
    ("N[C@H](Cc1ccc(F)c(F)c1)C(=O)O",         "D_Phe_34diF", "D-3,4-Difluorophenylalanine",           "aa","modified"),
    ("N[C@H](Cc1cc(F)cc(F)c1)C(=O)O",         "D_Phe_35diF", "D-3,5-Difluorophenylalanine",           "aa","modified"),
    ("N[C@H](Cc1c(F)cccc1F)C(=O)O",           "D_Phe_26diF", "D-2,6-Difluorophenylalanine",           "aa","modified"),
    ("N[C@H](Cc1ccc(F)cc1F)C(=O)O",           "D_Phe_24diF", "D-2,4-Difluorophenylalanine",           "aa","modified"),
    ("N[C@H](Cc1ccc(C(C)(C)C)cc1)C(=O)O",     "D_Phe_4tBu",  "D-4-tert-Butylphenylalanine",          "aa","modified"),
    ("N[C@H](Cc1c(F)c(F)c(F)c(F)c1F)C(=O)O", "D_Phe_F5",    "D-Pentafluorophenylalanine",             "aa","modified"),
    ("N[C@H](Cc1ccc(C(C)=O)cc1)C(=O)O",       "D_Phe_4Ac",   "D-4-Acetylphenylalanine",               "aa","modified"),
    ("N[C@H](Cc1ccccc1[N+](=O)[O-])C(=O)O",   "D_Phe_2NO2",  "D-2-Nitrophenylalanine",                "aa","modified"),
    ("N[C@H](Cc1cccc([N+](=O)[O-])c1)C(=O)O", "D_Phe_3NO2",  "D-3-Nitrophenylalanine",                "aa","modified"),
    ("N[C@H](Cc1cccc(O)c1)C(=O)O",            "D_Phe_3OH",   "D-3-Hydroxyphenylalanine",              "aa","modified"),

    # ── D-Trp ring variants ───────────────────────────────────────────────────
    ("N[C@H](Cc1c[nH]c2cc(F)ccc12)C(=O)O",    "D_Trp_5F",    "D-5-Fluorotryptophan",                  "aa","modified"),
    ("N[C@H](Cc1c[nH]c2cc(C)ccc12)C(=O)O",    "D_Trp_5Me",   "D-5-Methyltryptophan",                  "aa","modified"),
    ("N[C@H](Cc1c[nH]c2cc(Br)ccc12)C(=O)O",   "D_Trp_5Br",   "D-5-Bromotryptophan",                   "aa","modified"),
    ("N[C@H](Cc1c[nH]c2cc(OC)ccc12)C(=O)O",   "D_Trp_5OMe",  "D-5-Methoxytryptophan",                 "aa","modified"),
    ("N[C@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O",    "D_Trp_6F",    "D-6-Fluorotryptophan",                  "aa","modified"),
    ("N[C@H](Cc1c[nH]c2ncccc12)C(=O)O",       "D_Trp_7Aza",  "D-7-Azatryptophan",                     "aa","modified"),
    ("N[C@H](Cc1cn(C)c2ccccc12)C(=O)O",       "D_Trp_1Me",   "D-1-Methyltryptophan",                  "aa","modified"),
    ("N[C@H](Cc1c[nH]c2c(F)cccc12)C(=O)O",    "D_Trp_4F",    "D-4-Fluorotryptophan",                  "aa","modified"),

    # ── D-Tyr ring variants ───────────────────────────────────────────────────
    ("N[C@H](Cc1ccc(O)c(F)c1)C(=O)O",         "D_Tyr_3F",    "D-3-Fluorotyrosine",                    "aa","modified"),
    ("N[C@H](Cc1ccc(O)c(Cl)c1)C(=O)O",        "D_Tyr_3Cl",   "D-3-Chlorotyrosine",                    "aa","modified"),
    ("N[C@H](Cc1ccc(O)c(Br)c1)C(=O)O",        "D_Tyr_3Br",   "D-3-Bromotyrosine",                     "aa","modified"),
    ("N[C@H](Cc1c(C)cc(O)cc1C)C(=O)O",        "D_Tyr_35diMe","D-3,5-Dimethyltyrosine",                "aa","modified"),

    # ── D-Pro analogs ─────────────────────────────────────────────────────────
    ("OC(=O)[C@H]1C[C@@H](F)CN1",             "D_Pro_4F_t",  "D-trans-4-Fluoroproline",               "aa","modified"),
    ("OC(=O)[C@H]1C[C@H](F)CN1",              "D_Pro_4F_c",  "D-cis-4-Fluoroproline",                 "aa","modified"),
    ("OC(=O)[C@H]1Cc2ccccc2N1",               "D_Ind",       "D-Indoline-2-carboxylic acid",           "aa","modified"),

    # ── D-Glu/Asp side-chain esters ──────────────────────────────────────────
    ("N[C@H](CCC(=O)OC)C(=O)O",               "D_Glu_OMe",   "D-Glu gamma-methyl ester",              "aa","modified"),
    ("N[C@H](CCC(=O)OC(C)(C)C)C(=O)O",        "D_Glu_OtBu",  "D-Glu gamma-tert-butyl ester",         "aa","modified"),
    ("N[C@H](CCC(=O)OCc1ccccc1)C(=O)O",       "D_Glu_OBn",   "D-Glu gamma-benzyl ester",             "aa","modified"),
    ("N[C@H](CC(=O)OC(C)(C)C)C(=O)O",         "D_Asp_OtBu",  "D-Asp beta-tert-butyl ester",          "aa","modified"),
    ("N[C@H](CC(=O)OCc1ccccc1)C(=O)O",        "D_Asp_OBn",   "D-Asp beta-benzyl ester",              "aa","modified"),

    # ── More L-Trp ring variants ──────────────────────────────────────────────
    ("N[C@@H](Cc1c[nH]c2c(OC)cccc12)C(=O)O",  "Trp_4OMe",   "L-4-Methoxytryptophan",                 "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cc([N+](=O)[O-])ccc12)C(=O)O", "Trp_5NO2", "L-5-Nitrotryptophan",            "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2ccc(Cl)cc12)C(=O)O",  "Trp_6Cl",    "L-6-Chlorotryptophan",                  "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2ccc(Br)cc12)C(=O)O",  "Trp_6Br",    "L-6-Bromotryptophan",                   "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2ccc(C)cc12)C(=O)O",   "Trp_6Me",    "L-6-Methyltryptophan",                  "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O",   "Trp_6F2",    "L-6-Fluorotryptophan (alt regio)",      "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cc(Cl)ccc12)C(=O)O",  "Trp_5Cl",    "L-5-Chlorotryptophan",                  "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cc(I)ccc12)C(=O)O",   "Trp_5I",     "L-5-Iodotryptophan",                    "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2c([N+](=O)[O-])cccc12)C(=O)O", "Trp_4NO2", "L-4-Nitrotryptophan",            "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cccc(F)c12)C(=O)O",   "Trp_7F",     "L-7-Fluorotryptophan",                  "aa","modified"),

    # ── Unusual / natural-source AAs ─────────────────────────────────────────
    # Kynurenine: tryptophan ring-opening catabolite
    ("N[C@@H](CC(=O)c1ccccc1N)C(=O)O",        "Kyn",        "L-Kynurenine",                           "aa","modified"),

    # Homo-lysine (2,7-diaminoheptanoic acid)
    ("N[C@@H](CCCCCN)C(=O)O",                 "hLys",       "L-Homolysine",                           "aa","modified"),

    # Cyclobutylalanine / alpha-cyclobutyl-glycine
    ("N[C@@H](CC1CCC1)C(=O)O",                "cBuAla",     "L-Cyclobutylalanine",                    "aa","modified"),

    # Cyclopentylalanine
    ("N[C@@H](CC1CCCC1)C(=O)O",               "cPenAla",    "L-Cyclopentylalanine",                   "aa","modified"),

    # Furylalanine (2-furyl)
    ("N[C@@H](Cc1ccco1)C(=O)O",               "Fur2Ala",    "L-2-Furyl-alanine",                      "aa","modified"),

    # Furylalanine (3-furyl)
    ("N[C@@H](Cc1ccoc1)C(=O)O",               "Fur3Ala",    "L-3-Furyl-alanine",                      "aa","modified"),

    # Biphenylalanine (4-phenyl-L-Phe)
    ("N[C@@H](Cc1ccc(-c2ccccc2)cc1)C(=O)O",   "BiPheAla",   "L-4-Biphenylalanine",                   "aa","modified"),

    # Nε-Arginine(NO2) — protected Arg (nitroguanidine)
    ("N[C@@H](CCCNC(=N)N[N+](=O)[O-])C(=O)O", "ArgNO2",    "L-Arg(Nomega-nitro)",                    "aa","modified"),

    # trans-3-Hydroxyproline
    ("OC(=O)[C@@H]1C[C@H](O)CN1",             "Pro_3OH_t",  "trans-3-Hydroxyproline (2S,3R)",         "aa","modified"),

    # cis-3-Hydroxyproline
    ("OC(=O)[C@@H]1C[C@@H](O)CN1",            "Pro_3OH_c",  "cis-3-Hydroxyproline (2S,3S)",           "aa","modified"),

    # 3-Fluoroproline (racemate / unspecified stereo)
    ("OC(=O)C1CC(F)CN1",                      "Pro_3F",     "3-Fluoroproline",                        "aa","modified"),

    # Nε-formyl-L-lysine
    ("N[C@@H](CCCCNC=O)C(=O)O",              "Lys_For",    "N-epsilon-formyl-L-lysine",               "aa","modified"),

    # L-Biotin-Lys (Nε-biotinyl, streptavidin handle)
    ("N[C@@H](CCCCNC(=O)CCCCC1SCC2NC(=O)NC21)C(=O)O", "BiotinLys",
     "N-epsilon-biotinyl-L-lysine", "aa","modified"),

    # L-2-Amino-3-(4-azidophenyl)propanoic acid (azido-Phe, bioorthogonal)
    ("N[C@@H](Cc1ccc(N=[N+]=[N-])cc1)C(=O)O", "Phe_4Az",   "L-para-Azidophenylalanine",              "aa","modified"),

    # L-2-Amino-3-(4-propargyloxyphenyl)propanoic acid (O-propargyl-Tyr, alkyne handle)
    ("N[C@@H](Cc1ccc(OCC#C)cc1)C(=O)O",       "Phe_4OPrg", "L-para-O-Propargyl-Tyr",                "aa","modified"),

    # L-para-Cyanophenylalanine (already checked: Phe_4CN above, same as... wait that's D above)
    # (already have L-Phe_4CN from batch2) — skip

    # Dopa (3,4-dihydroxyphenylalanine) — same as L-DOPA
    ("N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O",        "Dopa",      "L-3,4-Dihydroxyphenylalanine (L-DOPA)",  "aa","modified"),

    # 4-Fluoromethyl-L-Phe (unusual CF2H/CH2F handle)
    ("N[C@@H](Cc1ccc(CF)cc1)C(=O)O",          "Phe_4CF",   "L-4-(Fluoromethyl)phenylalanine",        "aa","modified"),

    # ── More PEG linkers ─────────────────────────────────────────────────────
    # OEG2 (2-unit, 8-amino-3,6-dioxaoctanoic acid)
    ("NCCOCCOC(=O)O",                          "OEG2",      "8-Amino-3,6-dioxaoctanoic acid (mini-PEG 2-unit)", "linker","linker"),

    # OEG5 (5-unit)
    ("NCCOCCOCCOCCOCCOCC(=O)O",                "OEG5",      "17-Amino-3,6,9,12,15-pentaoxaheptadecanoic acid (mini-PEG 5-unit)", "linker","linker"),

    # OEG6 (6-unit)
    ("NCCOCCOCCOCCOCCOCCOCC(=O)O",             "OEG6",      "20-Amino-3,6,9,12,15,18-hexaoxaicosanoic acid (mini-PEG 6-unit)", "linker","linker"),

    # ── D-hPhe2 (D-homophenylalanine, 2-amino-4-phenylbutanoic) ─────────────
    ("N[C@H](CCc1ccccc1)C(=O)O",              "D_hPhe2",   "D-Homophenylalanine (2-amino-4-phenylbutanoic acid)", "aa","modified"),

    # ── More halogenated Phe ────────────────────────────────────────────────
    # L-Phe_4I (iodo already in library — checking: Phe_4I is in list, skip)
    # 3,4,5-Trifluoro-Phe (more symmetric than F5)
    ("N[C@@H](Cc1cc(F)c(F)c(F)c1)C(=O)O",    "Phe_345triF","L-3,4,5-Trifluorophenylalanine",         "aa","modified"),

    # ── D-β-amino acids ───────────────────────────────────────────────────────
    # D-forms of β³-homo AAs (swap [C@@H] to [C@H] in the CHUCKLES β-carbon)
    # These need manual CHUCKLES: [1*]N([3*])[C@H](sidechain)CC([2*])=O
    # Handled below as manual entries

    # ── Cyclic amino acids ────────────────────────────────────────────────────
    # 1-Aminocyclopropane-1-carboxylic acid (ACC)
    ("NCC1(CC1)C(=O)O",                       "ACC",       "1-Aminocyclopropane-1-carboxylic acid", "aa","modified"),

    # 1-Aminocyclobutane-1-carboxylic acid
    ("NC1(CCC1)C(=O)O",                       "ACBC",      "1-Aminocyclobutane-1-carboxylic acid", "aa","modified"),

    # 1-Aminocyclohexane-1-carboxylic acid (Ac6c is already there, skip)

    # ── N-alkyl glycines (peptoid building blocks) ────────────────────────────
    # Nspe (N-(S-1-phenylethyl)glycine — common peptoid)
    ("N([C@@H](C)c1ccccc1)CC(=O)O",           "Nspe",      "N-(S-1-phenylethyl)glycine (peptoid)", "aa","modified"),

    # Nrpe (N-(R-1-phenylethyl)glycine)
    ("N([C@H](C)c1ccccc1)CC(=O)O",            "Nrpe",      "N-(R-1-phenylethyl)glycine (peptoid)", "aa","modified"),

    # N-octylglycine (lipopeptoid handle)
    ("N(CCCCCCCC)CC(=O)O",                    "NocGly",    "N-octylglycine (lipopeptoid)",          "aa","modified"),

    # ── More aromatic / heterocyclic ─────────────────────────────────────────
    # Thiazolylalanine (2-thiazolyl)
    ("N[C@@H](Cc1nccs1)C(=O)O",               "Thz2Ala",   "L-2-(Thiazolyl)alanine",               "aa","modified"),

    # Pyrazol-3-yl-Ala
    ("N[C@@H](Cc1cc[nH]n1)C(=O)O",            "Pyz3Ala",   "L-3-(Pyrazolyl)alanine",               "aa","modified"),

    # Imidazol-2-yl-Ala (2-substituted His isomer)
    ("N[C@@H](Cc1ncc[nH]1)C(=O)O",            "Im2Ala",    "L-2-(Imidazolyl)alanine",              "aa","modified"),

    # ── Di-amino acids (branching points) ────────────────────────────────────
    # 2,3-Diaminopropionic acid (Dap already there ✓)
    # 2,4-Diaminobutyric acid (Dab already there ✓)
    # 2,5-Diaminopentanoic acid (Orn already there ✓)
    # 2,7-Diaminoheptanoic acid = hLys (added above)
    # 2,8-Diaminooctanoic acid
    ("N[C@@H](CCCCCCN)C(=O)O",                "Dao",       "L-2,8-Diaminooctanoic acid",           "aa","modified"),

    # ── Carbocyclic AA analogues ──────────────────────────────────────────────
    # ACPC (aminocyclopentane carboxylic acid, cis- scaffold)
    ("OC(=O)[C@@H]1CCC[C@H]1N",              "ACPC",      "cis-(1R,2S)-ACPC scaffold residue",    "aa","modified"),
]


# D-β³-amino acids: manual CHUCKLES only (pre_activate can't find β-backbone)
RGROUPS_STD = "[H],[OH],[H],None,None,None"
CT_STD = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

MANUAL_MONOMERS = [
    # (chuckles, abbr, name, m_type, m_subtype, rgroups, chem_types)
    ("[1*]N([3*])[C@H](Cc1ccccc1)CC([2*])=O",   "D_b3hPhe", "D-beta3-Homophenylalanine", "aa","modified", RGROUPS_STD, CT_STD),
    ("[1*]N([3*])[C@H](CC(C)C)CC([2*])=O",      "D_b3hLeu", "D-beta3-Homoleucine",       "aa","modified", RGROUPS_STD, CT_STD),
    ("[1*]N([3*])[C@H](C(C)C)CC([2*])=O",       "D_b3hVal", "D-beta3-Homovaline",        "aa","modified", RGROUPS_STD, CT_STD),
    ("[1*]N([3*])[C@H](Cc1c[nH]c2ccccc12)CC([2*])=O", "D_b3hTrp", "D-beta3-Homotryptophan", "aa","modified", RGROUPS_STD, CT_STD),
    ("[1*]N([3*])[C@H](CO)CC([2*])=O",          "D_b3hSer", "D-beta3-Homoserine",        "aa","modified", RGROUPS_STD, CT_STD),
    ("[1*]N([3*])[C@H](C)CC([2*])=O",           "D_b3hAla", "D-beta3-Homoalanine",       "aa","modified", RGROUPS_STD, CT_STD),
]


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor
    from pyPept.interfaces.monomer_pipeline import pre_activate

    rdDepictor.SetPreferCoordGen(True)

    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    existing = set()
    for m in suppl:
        if m:
            existing.add(m.GetPropsAsDict().get("m_abbr", ""))
    print(f"Existing monomers: {len(existing)}")

    added = skipped = 0
    failed = []

    def _fmt_rgroups(leaving_dict):
        slots = ["None"] * 6
        for slot, lg in leaving_dict.items():
            idx = int(slot) - 1
            if 0 <= idx < 6:
                slots[idx] = lg
        return ",".join(slots)

    def _fmt_cts(cts_dict):
        return ",".join(f"{k}:{v}" for k, v in sorted(cts_dict.items()))

    with open(SDF_PATH, "a", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        # pre_activate entries
        for smiles, abbr, name, mtype, msubtype in NEW_MONOMERS:
            if abbr in existing:
                skipped += 1
                continue
            try:
                result = pre_activate(smiles)
                mol = Chem.MolFromSmiles(result.chuckles)
                if mol is None:
                    raise ValueError("MolFromSmiles returned None on CHUCKLES")
                rdDepictor.Compute2DCoords(mol)
                mol.SetProp("m_abbr", abbr)
                mol.SetProp("symbol", abbr)
                mol.SetProp("m_name", name)
                mol.SetProp("m_type", mtype)
                mol.SetProp("m_subtype", msubtype)
                mol.SetProp("m_Rgroups", _fmt_rgroups(result.leaving))
                mol.SetProp("m_chem_types", _fmt_cts(result.chem_types))
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")
            except Exception as exc:
                failed.append((abbr, str(exc).split("\n")[0][:100]))

        # manual CHUCKLES entries
        for chuckles, abbr, name, mtype, msubtype, rgroups, chem_types in MANUAL_MONOMERS:
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
                mol.SetProp("m_Rgroups", rgroups)
                mol.SetProp("m_chem_types", chem_types)
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")
            except Exception as exc:
                failed.append((abbr, str(exc).split("\n")[0][:100]))

        writer.close()

    print(f"\nDone. Added {added} monomers (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")


if __name__ == "__main__":
    main()
