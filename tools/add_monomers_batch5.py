#!/usr/bin/env python3
"""Batch 5: ~200 new monomers.

Groups:
  - D-forms of batch4 unusual AAs
  - More Phe/Tyr/Trp ring variants (L and D)
  - Heteroaromatic AAs (BimAla, QuiAla, Thi3Ala, etc.)
  - More Pro ring variants
  - Homo-AA series (hLeu, hMet, hIle, hNle, hSer, hCys, hCit)
  - Misc useful AAs (allylglycine, cyanoalanine, azide handles, etc.)
  - RCM-staple AAs: S3, R3, R5, S8
  - NMe-unusual AAs (L and D)
  - Alpha-methyl AAs (missing from existing series)
  - D-alpha-methyl AAs
  - Orthogonal PG variants (Alloc-Lys, OAll-Glu, OAll-Asp, O-allyl Ser/Thr/Tyr/Cys)
  - Dehydro AAs (Dha, DeltaAbu, DeltaPhe)
  - Short-chain diacid caps (Succ, Glt, Adp, Pim, Sub, Azel, Seb)
  - Acyl caps (ClAc, BrAc, Propionyl, Butyryl)
  - Extra omega-amino-acid linkers (9, 10, 12 carbon)

Run from repo root:
    python tools/add_monomers_batch5.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# ── pre_activate entries (standard alpha-AAs) ────────────────────────────────
# Tuple: (smiles, abbr, name, m_type, m_subtype)
NEW_MONOMERS = [

    # ── D-forms of batch4 unusual AAs ────────────────────────────────────────
    ("N[C@H](Cc1ccco1)C(=O)O",                    "D_Fur2Ala",  "D-2-Furylalanine",                       "aa","modified"),
    ("N[C@H](Cc1ccoc1)C(=O)O",                    "D_Fur3Ala",  "D-3-Furylalanine",                       "aa","modified"),
    ("N[C@H](Cc1ccc(-c2ccccc2)cc1)C(=O)O",        "D_BiPheAla", "D-4-Biphenylalanine",                    "aa","modified"),
    ("N[C@H](CC1CCC1)C(=O)O",                     "D_cBuAla",   "D-Cyclobutylalanine",                    "aa","modified"),
    ("N[C@H](CC1CCCC1)C(=O)O",                    "D_cPenAla",  "D-Cyclopentylalanine",                   "aa","modified"),
    ("N[C@H](CC(=O)c1ccccc1N)C(=O)O",             "D_Kyn",      "D-Kynurenine",                           "aa","modified"),
    ("N[C@H](CCCNC(=N)N[N+](=O)[O-])C(=O)O",     "D_ArgNO2",   "D-Arg(Nomega-nitro)",                    "aa","modified"),
    ("N[C@H](Cc1ccc(O)c(O)c1)C(=O)O",             "D_Dopa",     "D-3,4-Dihydroxyphenylalanine (D-DOPA)",  "aa","modified"),
    ("N[C@H](CN=[N+]=[N-])C(=O)O",                "D_AzAla",    "D-Azidoalanine",                         "aa","modified"),
    ("N[C@H](CC1CC1)C(=O)O",                      "D_CyoAla",   "D-Cyclopropylalanine",                   "aa","modified"),
    ("N[C@H](Cc1nccs1)C(=O)O",                    "D_Thz2Ala",  "D-2-(Thiazolyl)alanine",                 "aa","modified"),
    ("N[C@H](Cc1cc[nH]n1)C(=O)O",                 "D_Pyz3Ala",  "D-3-(Pyrazolyl)alanine",                 "aa","modified"),
    ("N[C@H](Cc1ncc[nH]1)C(=O)O",                 "D_Im2Ala",   "D-2-(Imidazolyl)alanine",                "aa","modified"),
    ("N[C@H](CCCCNC(=O)CCCCC1SCC2NC(=O)NC21)C(=O)O", "D_BiotinLys", "N-epsilon-biotinyl-D-lysine",        "aa","modified"),
    ("N[C@H](CCCCCN)C(=O)O",                      "D_hLys",     "D-Homolysine",                           "aa","modified"),
    ("N[C@H](CCCCCCN)C(=O)O",                     "D_Dao",      "D-2,8-Diaminooctanoic acid",             "aa","modified"),
    ("N[C@H](Cc1ccc(N=[N+]=[N-])cc1)C(=O)O",      "D_Phe_4Az",  "D-para-Azidophenylalanine",              "aa","modified"),

    # ── More L-Phe ring variants ──────────────────────────────────────────────
    ("N[C@@H](Cc1ccc2c(c1)OCO2)C(=O)O",           "PipAla",     "L-Piperonylalanine (3,4-methylenedioxy-Phe)", "aa","modified"),
    ("N[C@@H](Cc1ccc(N(C)C)cc1)C(=O)O",           "Phe_4NMe2",  "L-4-(Dimethylamino)phenylalanine",       "aa","modified"),
    ("N[C@@H](Cc1ccc(C(C)C)cc1)C(=O)O",           "Phe_4iPr",   "L-4-Isopropylphenylalanine",             "aa","modified"),
    ("N[C@@H](Cc1ccc(OC)c(OC)c1)C(=O)O",          "Phe_34diOMe","L-3,4-Dimethoxyphenylalanine",           "aa","modified"),
    ("N[C@@H](Cc1ccc(N)cc1)C(=O)O",               "Phe_4NH2",   "L-4-Aminophenylalanine",                 "aa","modified"),
    ("N[C@@H](Cc1ccc(Br)cc1)C(=O)O",              "Phe_4Br",    "L-4-Bromophenylalanine",                 "aa","modified"),
    ("N[C@@H](Cc1ccc(I)cc1)C(=O)O",               "Phe_4I",     "L-4-Iodophenylalanine",                  "aa","modified"),
    ("N[C@@H](Cc1ccc(SC)cc1)C(=O)O",              "Phe_4SMe",   "L-4-(Methylthio)phenylalanine",          "aa","modified"),
    ("N[C@@H](Cc1ccc([N+](=O)[O-])cc1)C(=O)O",    "Phe_4NO2",   "L-4-Nitrophenylalanine",                 "aa","modified"),

    # ── D-Phe more variants ───────────────────────────────────────────────────
    ("N[C@H](Cc1ccc2c(c1)OCO2)C(=O)O",            "D_PipAla",   "D-Piperonylalanine",                     "aa","modified"),
    ("N[C@H](Cc1ccc(Br)cc1)C(=O)O",               "D_Phe_4Br",  "D-4-Bromophenylalanine",                 "aa","modified"),
    ("N[C@H](Cc1ccc(I)cc1)C(=O)O",                "D_Phe_4I",   "D-4-Iodophenylalanine",                  "aa","modified"),
    ("N[C@H](Cc1ccc([N+](=O)[O-])cc1)C(=O)O",     "D_Phe_4NO2", "D-4-Nitrophenylalanine",                 "aa","modified"),
    ("N[C@H](Cc1ccc(N)cc1)C(=O)O",                "D_Phe_4NH2", "D-4-Aminophenylalanine",                 "aa","modified"),
    ("N[C@H](Cc1ccc(N(C)C)cc1)C(=O)O",            "D_Phe_4NMe2","D-4-(Dimethylamino)phenylalanine",       "aa","modified"),
    ("N[C@H](Cc1ccc(SC)cc1)C(=O)O",               "D_Phe_4SMe", "D-4-(Methylthio)phenylalanine",          "aa","modified"),

    # ── Heteroaromatic AAs (L) ────────────────────────────────────────────────
    ("N[C@@H](Cc1nc2ccccc2[nH]1)C(=O)O",          "BimAla",     "L-Benzimidazol-2-yl-alanine",            "aa","modified"),
    ("N[C@@H](Cc1ccc2ccccc2n1)C(=O)O",            "QuiAla",     "L-Quinolin-2-yl-alanine",                "aa","modified"),
    ("N[C@@H](Cc1nccc2ccccc12)C(=O)O",            "IsoQuiAla",  "L-Isoquinolin-1-yl-alanine",             "aa","modified"),
    ("N[C@@H](Cc1ccsc1)C(=O)O",                   "Thi3Ala",    "L-Thiophen-3-yl-alanine",                "aa","modified"),
    ("N[C@@H](Cc1[nH]nc2ccccc12)C(=O)O",          "IndAzAla",   "L-Indazol-3-yl-alanine",                 "aa","modified"),
    ("N[C@@H](Cn1cccc1)C(=O)O",                   "Pyr1Ala",    "L-Pyrrol-1-yl-alanine (N-linked)",       "aa","modified"),
    ("N[C@@H](Cn1ccnn1)C(=O)O",                   "TriazAla",   "L-1,2,3-Triazol-1-yl-alanine",           "aa","modified"),
    ("N[C@@H](Cc1nc2ccccc2o1)C(=O)O",             "BenzoxAla",  "L-Benzoxazol-2-yl-alanine",              "aa","modified"),

    # ── Heteroaromatic AAs (D) ────────────────────────────────────────────────
    ("N[C@H](Cc1nc2ccccc2[nH]1)C(=O)O",           "D_BimAla",   "D-Benzimidazol-2-yl-alanine",            "aa","modified"),
    ("N[C@H](Cc1ccc2ccccc2n1)C(=O)O",             "D_QuiAla",   "D-Quinolin-2-yl-alanine",                "aa","modified"),
    ("N[C@H](Cc1ccsc1)C(=O)O",                    "D_Thi3Ala",  "D-Thiophen-3-yl-alanine",                "aa","modified"),
    ("N[C@H](Cc1[nH]nc2ccccc12)C(=O)O",           "D_IndAzAla", "D-Indazol-3-yl-alanine",                 "aa","modified"),

    # ── More L-Trp ring variants ──────────────────────────────────────────────
    ("N[C@@H](Cc1c[nH]c2c(C)cccc12)C(=O)O",       "Trp_4Me",    "L-4-Methyltryptophan",                   "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2c(Cl)cccc12)C(=O)O",      "Trp_4Cl",    "L-4-Chlorotryptophan",                   "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2c(Br)cccc12)C(=O)O",      "Trp_4Br",    "L-4-Bromotryptophan",                    "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cc(N)ccc12)C(=O)O",       "Trp_5NH2",   "L-5-Aminotryptophan",                    "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cccc(C)c12)C(=O)O",       "Trp_7Me",    "L-7-Methyltryptophan",                   "aa","modified"),
    ("N[C@@H](Cc1c[nH]c2cccc(Cl)c12)C(=O)O",      "Trp_7Cl",    "L-7-Chlorotryptophan",                   "aa","modified"),
    ("N[C@@H](Cc1c(C)[nH]c2ccccc12)C(=O)O",       "Trp_2Me",    "L-2-Methyltryptophan",                   "aa","modified"),

    # ── D-Trp ring variants (additional) ─────────────────────────────────────
    ("N[C@H](Cc1c[nH]c2c(C)cccc12)C(=O)O",        "D_Trp_4Me",  "D-4-Methyltryptophan",                   "aa","modified"),
    ("N[C@H](Cc1c[nH]c2c(Br)cccc12)C(=O)O",       "D_Trp_4Br",  "D-4-Bromotryptophan",                    "aa","modified"),
    ("N[C@H](Cc1c[nH]c2cccc(C)c12)C(=O)O",        "D_Trp_7Me",  "D-7-Methyltryptophan",                   "aa","modified"),

    # ── Tyr variants ──────────────────────────────────────────────────────────
    ("N[C@@H](Cc1ccc(O)c(I)c1)C(=O)O",            "Tyr_3I",     "L-3-Iodotyrosine",                       "aa","modified"),
    ("N[C@@H](Cc1cc(I)c(O)c(I)c1)C(=O)O",         "Tyr_35diI",  "L-3,5-Diiodotyrosine",                   "aa","modified"),
    ("N[C@@H](Cc1ccc(O)c([N+](=O)[O-])c1)C(=O)O", "Tyr_3NO2",   "L-3-Nitrotyrosine",                      "aa","modified"),
    ("N[C@H](Cc1ccc(O)c(I)c1)C(=O)O",             "D_Tyr_3I",   "D-3-Iodotyrosine",                       "aa","modified"),
    ("N[C@H](Cc1ccc(O)c([N+](=O)[O-])c1)C(=O)O",  "D_Tyr_3NO2", "D-3-Nitrotyrosine",                      "aa","modified"),

    # ── Protected Tyr/Ser ─────────────────────────────────────────────────────
    ("N[C@@H](Cc1ccc(OCc2ccccc2)cc1)C(=O)O",      "Tyr_OBn",    "L-O-Benzyl-tyrosine",                    "aa","modified"),
    ("N[C@@H](COCc1ccccc1)C(=O)O",                "Ser_OBn",    "L-O-Benzyl-serine",                      "aa","modified"),

    # ── Pro analogs ───────────────────────────────────────────────────────────
    ("OC(=O)[C@@H]1C[C@H](C)CN1",                 "Pro_4Me_t",  "trans-4-Methyl-L-proline (2S,4R)",       "aa","modified"),
    ("OC(=O)[C@@H]1C[C@@H](C)CN1",                "Pro_4Me_c",  "cis-4-Methyl-L-proline (2S,4S)",         "aa","modified"),
    ("OC(=O)[C@@H]1C[C@H](N)CN1",                 "Pro_4NH2_t", "trans-4-Amino-L-proline (2S,4R)",        "aa","modified"),
    ("OC(=O)[C@@H]1C[C@@H](N)CN1",                "Pro_4NH2_c", "cis-4-Amino-L-proline (2S,4S)",          "aa","modified"),

    # ── Homo-amino acids (sidechain-extended) ────────────────────────────────
    ("N[C@@H](CCC(C)C)C(=O)O",                    "hLeu",       "L-Homoleucine (2-amino-5-methylhexanoic acid)", "aa","modified"),
    ("N[C@@H](CCCSC)C(=O)O",                      "hMet",       "L-Homomethionine (2-amino-5-methylthiopentanoic acid)", "aa","modified"),
    ("N[C@@H](CC(C)CC)C(=O)O",                    "hIle",       "L-Homoisoleucine (2-amino-3-methylhexanoic acid)", "aa","modified"),
    ("N[C@@H](CCCCC)C(=O)O",                      "hNle",       "L-Homon­orleucine (2-aminoheptanoic acid)", "aa","modified"),
    ("N[C@@H](CCO)C(=O)O",                        "hSer",       "L-Homoserine (2-amino-4-hydroxybutanoic acid)", "aa","modified"),
    ("N[C@@H](CCS)C(=O)O",                        "hCys",       "L-Homocysteine",                         "aa","modified"),
    ("N[C@@H](CCCCNC(=O)N)C(=O)O",               "hCit",       "L-Homocitrulline (Nomega-ureido)",         "aa","modified"),

    # ── D-Homo-amino acids ────────────────────────────────────────────────────
    ("N[C@H](CCC(C)C)C(=O)O",                     "D_hLeu",     "D-Homoleucine",                          "aa","modified"),
    ("N[C@H](CCCSC)C(=O)O",                       "D_hMet",     "D-Homomethionine",                       "aa","modified"),
    ("N[C@H](CC(C)CC)C(=O)O",                     "D_hIle",     "D-Homoisoleucine",                       "aa","modified"),
    ("N[C@H](CCCCC)C(=O)O",                       "D_hNle",     "D-Homonorleucine",                       "aa","modified"),
    ("N[C@H](CCO)C(=O)O",                         "D_hSer",     "D-Homoserine",                           "aa","modified"),
    ("N[C@H](CCS)C(=O)O",                         "D_hCys",     "D-Homocysteine",                         "aa","modified"),

    # ── Misc useful AAs ───────────────────────────────────────────────────────
    ("N[C@@H](CC=C)C(=O)O",                       "AllGly",     "L-Allylglycine (2-amino-4-pentenoic acid)", "aa","modified"),
    ("N[C@H](CC=C)C(=O)O",                        "D_AllGly",   "D-Allylglycine",                         "aa","modified"),
    ("N[C@@H](CC#N)C(=O)O",                       "CyanoAla",   "L-3-Cyanoalanine (beta-cyanoalanine)",   "aa","modified"),
    ("N[C@@H](CCCN=[N+]=[N-])C(=O)O",             "AzOrn",      "L-delta-Azidoornithine",                 "aa","modified"),
    ("N[C@H](CCCN=[N+]=[N-])C(=O)O",              "D_AzOrn",    "D-delta-Azidoornithine",                 "aa","modified"),
    ("N[C@@H](CCCCN=[N+]=[N-])C(=O)O",            "AzNva",      "L-epsilon-Azido-norvaline (5-azido)",    "aa","modified"),
    ("N[C@H](CCCCN=[N+]=[N-])C(=O)O",             "D_AzNva",    "D-epsilon-Azido-norvaline",              "aa","modified"),
    ("N[C@@H](CCCNC(=N)N)C(=O)O",                 "Nrg",        "L-Norarginine (4-guanidinobutyric AA)",  "aa","modified"),
    ("N[C@@H](CCCNC(=N)N(C)C)C(=O)O",             "ADMA",       "L-ADMA (asymmetric dimethylarginine)",   "aa","modified"),
    ("N[C@@H](CCCCNC)C(=O)O",                     "Lys_5Me",    "N-epsilon-methyl-L-lysine",              "aa","modified"),
    ("N[C@@H](CCCCN(C)C)C(=O)O",                  "Lys_5diMe",  "N-epsilon,epsilon-dimethyl-L-lysine",    "aa","modified"),
    ("N[C@@H](CCCCNC(=O)Oc1ccccc1)C(=O)O",        "Lys_Cbz",    "N-epsilon-Cbz-L-lysine",                "aa","modified"),
    ("N[C@@H](CCCNC(=O)Oc1ccccc1)C(=O)O",         "Orn_Cbz",    "N-delta-Cbz-L-ornithine",               "aa","modified"),
]


# ── Manual CHUCKLES entries ───────────────────────────────────────────────────
# Tuple: (chuckles, abbr, name, m_type, m_subtype, m_Rgroups, m_chem_types)

_AA_RG    = "[H],[OH],[H],None,None,None"   # standard AA (R1, R2, R3)
_NMEAA_RG = "[H],[OH],None,None,None,None"  # NMe-AA (R1, R2 only)
_STPL_RG  = "[H],[OH],[H],[H],None,None"    # staple AA (R1-R4)
_CAP_RG   = "None,[OH],None,None,None,None" # cap (R2 only)
_CAP4_RG  = "None,[OH],None,None,None,None" # same as above
_CT_AA    = "1:backbone_n,2:backbone_c,3:backbone_n_mod"
_CT_NME   = "1:backbone_n,2:backbone_c"
_CT_CAP   = ""
_CT_CAPC  = "2:backbone_c"
_CT_LINK  = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

MANUAL_MONOMERS = [

    # ── RCM staple AAs (Verdine hydrocarbon staple series) ───────────────────
    # [4*] marks the terminal alkene carbon; [3*] = N-mod; [@@]/[C@] = (S)/(R)
    ("[1*]N([3*])[C@@](C)(CC([4*])=C)C([2*])=O",          "S3", "(S)-alpha-Me-allylglycine (i,i+3 RCM)",  "aa","modified", _STPL_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CC([4*])=C)C([2*])=O",           "R3", "(R)-alpha-Me-allylglycine (i,i+3 RCM)",  "aa","modified", _STPL_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CCCC([4*])=C)C([2*])=O",         "R5", "(R)-alpha-Me-pentenylglycine (i,i+4 RCM, S5 mirror)", "aa","modified", _STPL_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(CCCCCCC([4*])=C)C([2*])=O",     "S8", "(S)-alpha-Me-octenylglycine (i,i+7 RCM, R8 mirror)", "aa","modified", _STPL_RG, _CT_AA),

    # ── NMe-unusual L-AAs ─────────────────────────────────────────────────────
    ("[1*]N(C)[C@@H](CC1CCCCC1)C([2*])=O",                "meCha",      "NMe-L-cyclohexylalanine",           "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](CCCC)C([2*])=O",                     "meNle",      "NMe-L-norleucine",                  "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](CCC)C([2*])=O",                      "meNva",      "NMe-L-norvaline",                   "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](CC)C([2*])=O",                       "meAbu",      "NMe-L-Abu (NMe-alpha-aminobutyric)", "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1cccc2ccccc12)C([2*])=O",          "me1Nal",     "NMe-L-1-naphthylalanine",           "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc2ccccc2c1)C([2*])=O",          "me2Nal",     "NMe-L-2-naphthylalanine",           "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(F)cc1)C([2*])=O",             "mePhe_4F",   "NMe-L-4-fluorophenylalanine",       "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(Cl)cc1)C([2*])=O",            "mePhe_4Cl",  "NMe-L-4-chlorophenylalanine",       "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(OC)cc1)C([2*])=O",            "mePhe_4OMe", "NMe-L-4-methoxyphenylalanine",      "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(Br)cc1)C([2*])=O",            "mePhe_4Br",  "NMe-L-4-bromophenylalanine",        "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1cccc(F)c1)C([2*])=O",             "mePhe_3F",   "NMe-L-3-fluorophenylalanine",       "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1cccc(OC)c1)C([2*])=O",            "mePhe_3OMe", "NMe-L-3-methoxyphenylalanine",      "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(-c2ccccc2)cc1)C([2*])=O",     "meBip",      "NMe-L-biphenylalanine",             "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](CC(c1ccccc1)c1ccccc1)C([2*])=O",     "meDip",      "NMe-L-3,3-diphenylalanine",         "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](C1CCCCC1)C([2*])=O",                 "meChg",      "NMe-L-cyclohexylglycine",           "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccncc1)C([2*])=O",                "me4Pal",     "NMe-L-4-pyridylalanine",            "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1cccnc1)C([2*])=O",                "me3Pal",     "NMe-L-3-pyridylalanine",            "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1c[nH]c2cc(F)ccc12)C([2*])=O",    "meTrp_5F",   "NMe-L-5-fluorotryptophan",          "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1c[nH]c2cc(C)ccc12)C([2*])=O",    "meTrp_5Me",  "NMe-L-5-methyltryptophan",          "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1c[nH]c2ncccc12)C([2*])=O",       "meTrp_7Aza", "NMe-L-7-azatryptophan",             "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(C(F)(F)F)cc1)C([2*])=O",     "mePhe_4CF3", "NMe-L-4-(trifluoromethyl)phenylalanine", "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@@H](Cc1ccc(C(C)(C)C)cc1)C([2*])=O",     "mePhe_4tBu", "NMe-L-4-tert-butylphenylalanine",   "aa","modified", _NMEAA_RG, _CT_NME),

    # ── NMe-D-AAs ─────────────────────────────────────────────────────────────
    ("[1*]N(C)[C@H](C)C([2*])=O",                         "D_meA",  "NMe-D-alanine",                     "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](Cc1ccccc1)C([2*])=O",                 "D_meF",  "NMe-D-phenylalanine",               "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](CC(C)C)C([2*])=O",                    "D_meL",  "NMe-D-leucine",                     "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](C(C)C)C([2*])=O",                     "D_meV",  "NMe-D-valine",                      "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H]([C@@H](C)CC)C([2*])=O",               "D_meI",  "NMe-D-isoleucine",                  "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](Cc1c[nH]c2ccccc12)C([2*])=O",         "D_meW",  "NMe-D-tryptophan",                  "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](Cc1ccc(O)cc1)C([2*])=O",              "D_meY",  "NMe-D-tyrosine",                    "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](CO)C([2*])=O",                        "D_meS",  "NMe-D-serine",                      "aa","modified", _NMEAA_RG, _CT_NME),
    ("[1*]N(C)[C@H](CCCC)C([2*])=O",                      "D_meNle",    "NMe-D-norleucine",              "aa","modified", _NMEAA_RG, _CT_NME),

    # ── Alpha-methyl AAs (missing from existing series) ───────────────────────
    ("[1*]N([3*])[C@@](C)(CCSC)C([2*])=O",                "aMeMet",  "L-alpha-methylmethionine (S)",      "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(CCCC)C([2*])=O",                "aMeNle",  "L-alpha-methylnorleucine (S)",      "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(CCC)C([2*])=O",                 "aMeNva",  "L-alpha-methylnorvaline (S)",       "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(CC)C([2*])=O",                  "aMeAbu",  "L-alpha-methyl-Abu (S)",            "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)([C@@H](C)O)C([2*])=O",          "aMeThr",  "L-alpha-methylthreonine (S)",       "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(Cc1ccc(F)cc1)C([2*])=O",        "aMePhe_4F",  "L-alpha-methyl-4-F-Phe (S)",    "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(Cc1ccc(Cl)cc1)C([2*])=O",       "aMePhe_4Cl", "L-alpha-methyl-4-Cl-Phe (S)",   "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@](C)(Cc1ccc(OC)cc1)C([2*])=O",       "aMePhe_4OMe","L-alpha-methyl-4-OMe-Phe (S)",  "aa","modified", _AA_RG, _CT_AA),

    # ── D-alpha-methyl AAs ────────────────────────────────────────────────────
    ("[1*]N([3*])[C@](C)(Cc1ccccc1)C([2*])=O",            "D_aMePhe",  "D-alpha-methylphenylalanine (R)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CC(C)C)C([2*])=O",               "D_aMeLeu",  "D-alpha-methylleucine (R)",       "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(C(C)C)C([2*])=O",                "D_aMeVal",  "D-alpha-methylvaline (R)",        "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)([C@@H](C)CC)C([2*])=O",          "D_aMeIle",  "D-alpha-methylisoleucine (R)",    "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(Cc1c[nH]c2ccccc12)C([2*])=O",    "D_aMeTrp",  "D-alpha-methyltryptophan (R)",    "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CCC(=O)O)C([2*])=O",             "D_aMeGlu",  "D-alpha-methylglutamic acid (R)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CC(=O)O)C([2*])=O",              "D_aMeAsp",  "D-alpha-methylaspartic acid (R)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(Cc1ccc(O)cc1)C([2*])=O",         "D_aMeTyr",  "D-alpha-methyltyrosine (R)",      "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CCC(=O)N)C([2*])=O",             "D_aMeGln",  "D-alpha-methylglutamine (R)",     "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CC(=O)N)C([2*])=O",              "D_aMeAsn",  "D-alpha-methylasparagine (R)",    "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@](C)(CO)C([2*])=O",                   "D_aMeSer",  "D-alpha-methylserine (R)",        "aa","modified", _AA_RG, _CT_AA),

    # ── Orthogonal protecting group variants ─────────────────────────────────
    # Alloc on amines (allyloxycarbonyl: NHC(=O)OCC=C)
    ("[1*]N([3*])[C@@H](CCCCNC(=O)OCC=C)C([2*])=O",      "Lys_Alloc",  "N-epsilon-Alloc-L-lysine",         "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](CCCNC(=O)OCC=C)C([2*])=O",       "Orn_Alloc",  "N-delta-Alloc-L-ornithine",        "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](CCNC(=O)OCC=C)C([2*])=O",        "Dab_Alloc",  "N-gamma-Alloc-L-Dab",              "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](CNC(=O)OCC=C)C([2*])=O",         "Dap_Alloc",  "N-beta-Alloc-L-Dap",               "aa","modified", _AA_RG, _CT_AA),
    # Allyl esters on acidic sidechains
    ("[1*]N([3*])[C@@H](CCC(=O)OCC=C)C([2*])=O",         "Glu_OAll",   "L-Glu gamma-allyl ester",          "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](CC(=O)OCC=C)C([2*])=O",          "Asp_OAll",   "L-Asp beta-allyl ester",           "aa","modified", _AA_RG, _CT_AA),
    # O-allyl ethers on hydroxyls
    ("[1*]N([3*])[C@@H](COCC=C)C([2*])=O",               "Ser_All",    "L-O-allylserine",                  "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H]([C@@H](C)OCC=C)C([2*])=O",       "Thr_All",    "L-O-allylthreonine",               "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](Cc1ccc(OCC=C)cc1)C([2*])=O",     "Tyr_All",    "L-O-allyltyrosine",                "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](CSCC=C)C([2*])=O",               "Cys_All",    "L-S-allylcysteine",                "aa","modified", _AA_RG, _CT_AA),
    # Mtt on lysine
    ("[1*]N([3*])[C@@H](CCCCNC(c1ccccc1)(c1ccccc1)c1ccc(C)cc1)C([2*])=O", "Lys_Mtt", "N-epsilon-Mtt-L-lysine", "aa","modified", _AA_RG, _CT_AA),

    # ── Dehydro AAs ───────────────────────────────────────────────────────────
    ("[1*]N([3*])C(=C)C([2*])=O",                         "Dha",       "Dehydroalanine (alpha,beta-didehydroalanine)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])C(=CC)C([2*])=O",                        "DeltaAbu",  "alpha,beta-Didehydroaminobutyric acid",        "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])C(=Cc1ccccc1)C([2*])=O",                 "DeltaPhe",  "alpha,beta-Didehydrophenylalanine",            "aa","modified", _AA_RG, _CT_AA),

    # ── β3-homo Phe/Trp ring variants (L and D) ───────────────────────────────
    ("[1*]N([3*])[C@@H](Cc1ccc(F)cc1)CC([2*])=O",         "b3hPhe_4F",  "beta3-Homo-4-fluorophenylalanine (L)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](Cc1ccc(Cl)cc1)CC([2*])=O",        "b3hPhe_4Cl", "beta3-Homo-4-chlorophenylalanine (L)", "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](Cc1ccc(OC)cc1)CC([2*])=O",        "b3hPhe_4OMe","beta3-Homo-4-methoxyphenylalanine (L)","aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](Cc1ccc(Br)cc1)CC([2*])=O",        "b3hPhe_4Br", "beta3-Homo-4-bromophenylalanine (L)",  "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@@H](Cc1c[nH]c2cc(F)ccc12)CC([2*])=O", "b3hTrp_5F",  "beta3-Homo-5-fluorotryptophan (L)",   "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@H](Cc1ccc(F)cc1)CC([2*])=O",          "D_b3hPhe_4F",  "D-beta3-Homo-4-fluorophenylalanine",  "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@H](Cc1ccc(Cl)cc1)CC([2*])=O",         "D_b3hPhe_4Cl", "D-beta3-Homo-4-chlorophenylalanine",  "aa","modified", _AA_RG, _CT_AA),
    ("[1*]N([3*])[C@H](Cc1ccc(OC)cc1)CC([2*])=O",         "D_b3hPhe_4OMe","D-beta3-Homo-4-methoxyphenylalanine", "aa","modified", _AA_RG, _CT_AA),
]


# ── Cap monomers (manual) ─────────────────────────────────────────────────────
# Tuple: (chuckles, abbr, name, m_Rgroups, m_chem_types)

CAP_MONOMERS = [
    # Diacid caps (short-chain, C-terminal ester attachment)
    ("[2*]OC(=O)CCC(=O)O",                "Succ",     "Succinyl cap (C4 diacid)",   _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCC(=O)O",               "Glt",      "Glutaryl cap (C5 diacid)",   _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCC(=O)O",              "Adp",      "Adipoyl cap (C6 diacid)",    _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCCC(=O)O",             "Pim",      "Pimeloyl cap (C7 diacid)",   _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCCCC(=O)O",            "Sub",      "Suberoyl cap (C8 diacid)",   _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCCCCC(=O)O",           "Azel",     "Azelaoyl cap (C9 diacid)",   _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCCCCCC(=O)O",          "Seb",      "Sebacoyl cap (C10 diacid)",  _CAP_RG, _CT_CAPC),
    ("[2*]OC(=O)CCCCCCCCCCC(=O)O",        "Dod_dAcid","Dodecanedioyl cap (C12 diacid)", _CAP_RG, _CT_CAPC),
    # Acyl caps (N-terminal)
    ("[2*]C(=O)CCl",                      "ClAc",     "Chloroacetyl cap",            _CAP_RG, _CT_CAP),
    ("[2*]C(=O)CBr",                      "BrAc",     "Bromoacetyl cap",             _CAP_RG, _CT_CAP),
    ("[2*]C(=O)CC",                       "PropCap",  "Propionyl cap",               _CAP_RG, _CT_CAP),
    ("[2*]C(=O)CCC",                      "ButCap",   "Butyryl cap",                 _CAP_RG, _CT_CAP),
    ("[2*]C(=O)CCCCC",                    "HexCap",   "Hexanoyl cap",                _CAP_RG, _CT_CAP),
    ("[2*]C(=O)C(F)(F)F",                 "CF3Ac",    "Trifluoroacetyl cap",         _CAP_RG, _CT_CAP),
]


# ── Linker monomers ───────────────────────────────────────────────────────────
# Tuple: (chuckles, abbr, name)  — type/subtype always linker/linker

LINKER_MONOMERS = [
    ("[1*]N([3*])CCCCCCCCC([2*])=O",      "9Anac",   "9-Aminononanoic acid linker"),
    ("[1*]N([3*])CCCCCCCCCC([2*])=O",     "10Adec",  "10-Aminodecanoic acid linker"),
    ("[1*]N([3*])CCCCCCCCCCCC([2*])=O",   "12Adod",  "12-Aminododecanoic acid linker"),
    ("[1*]N([3*])CCOCCOCCOCCOCCOCCOCCOCC([2*])=O", "OEG7", "7-unit PEG amino acid linker"),
    ("[1*]N([3*])CCOCCOCCOCCOCCOCCOCCOCCOC C([2*])=O", "OEG8_BAD", "placeholder-skip"),  # space intentional to trigger skip
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

        # ── pre_activate entries ──────────────────────────────────────────────
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
                failed.append((abbr, str(exc).split("\n")[0][:120]))

        # ── manual CHUCKLES entries ───────────────────────────────────────────
        for chuckles, abbr, name, mtype, msubtype, rgroups, chem_types in MANUAL_MONOMERS:
            if abbr.startswith("_skip") or abbr in existing:
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
                failed.append((abbr, str(exc).split("\n")[0][:120]))

        # ── cap monomers ──────────────────────────────────────────────────────
        for chuckles, abbr, name, rgroups, chem_types in CAP_MONOMERS:
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
                mol.SetProp("m_type", "cap")
                mol.SetProp("m_subtype", "cap")
                mol.SetProp("m_Rgroups", rgroups)
                mol.SetProp("m_chem_types", chem_types)
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")
            except Exception as exc:
                failed.append((abbr, str(exc).split("\n")[0][:120]))

        # ── linker monomers ───────────────────────────────────────────────────
        for chuckles, abbr, name in LINKER_MONOMERS:
            if "BAD" in abbr or abbr in existing:
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
                mol.SetProp("m_type", "linker")
                mol.SetProp("m_subtype", "linker")
                mol.SetProp("m_Rgroups", _AA_RG)
                mol.SetProp("m_chem_types", _CT_LINK)
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")
            except Exception as exc:
                failed.append((abbr, str(exc).split("\n")[0][:120]))

        writer.close()

    print(f"\nDone. Added {added} monomers (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")


if __name__ == "__main__":
    main()
