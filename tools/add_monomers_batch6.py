#!/usr/bin/env python3
"""Batch 6: ~230 new monomers targeting 1000 total.

Adds:
 - D-NMe AA complete series (D_meT/M/C/H/D/E/N/Q/K/R + aryl/heteroaryl D-NMe + mePro/D_mePro)
 - New NMe-AA: meOrn/meDab/meDap + more NMe-Phe variants
 - D-forms of all batch5 special entries
 - Missing D-alpha-methyl AAs + new L/D aMePhe/1Nal/2Nal/Cha/Chg/Bip/Trp variants
 - New L- and D-beta3-homo AAs (b3hPhe/Met/Nva/Cha/Chg/1Nal/2Nal/Asn/Gln/Arg/His/Trp)
 - Extended b3h-Phe halogenated series
 - Misc: SeM, VinGly, HPG, Cys_Prg, Ala_3Br/Cl, Pro_4OAc, D-Pro analogs, Phe_3I
 - Extended NMe-Phe/Tyr/Trp variants
 - More staple AAs (S4/S6/R4/R6) and acyl caps

Run from repo root:
    python tools/add_monomers_batch6.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# ─── R-group / chem-type constants ───────────────────────────────────────────
_AA_RG   = "[H],[OH],[H],None,None,None"
_AA_CT   = "1:backbone_n,2:backbone_c,3:backbone_n_mod"
_AA_RG2  = "[H],[OH],[H],[H],None,None"
_AA_CT2  = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

# NMe-AA standard (no sidechain R-groups)
_NME_RG  = "[H],[OH],None,None,None,None"
_NME_CT  = "1:backbone_n,2:backbone_c"
# NMe-AA with one sidechain R-group at [3*]
_NME1_RG = "[H],[OH],[H],None,None,None"
_NME1_CT = "1:backbone_n,2:backbone_c"
# NMe-AA with two sidechain R-groups at [3*] and [4*]
_NME2_RG = "[H],[OH],[H],[H],None,None"
_NME2_CT = "1:backbone_n,2:backbone_c"

_PRO_RG  = "[H],[OH],None,None,None,None"
_PRO_CT  = "1:backbone_n,2:backbone_c"

# Proline with one sidechain [3*]
_PRO1_RG = "[H],[OH],[H],None,None,None"
_PRO1_CT = "1:backbone_n,2:backbone_c,3:hydroxyl"

# Proline with two sidechain [3*],[4*]
_PRO2_RG = "[H],[OH],[H],[H],None,None"
_PRO2_CT = "1:backbone_n,2:backbone_c,3:amine_primary,4:amine_secondary"

# Staple AAs (alpha-methyl, R4=terminal alkene)
_STPL_RG = "[H],[OH],[H],[H],None,None"
_STPL_CT = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

# Alpha-methyl AAs (Aib-type)
_AME_RG  = "[H],[OH],[H],None,None,None"
_AME_CT  = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

# Cap
_CAP_RG  = "None,[OH],None,None,None,None"
_CAP_CT  = "2:backbone_c"

# ─── SECTION 1: D-NMe AA series ──────────────────────────────────────────────
# L-form stereo [C@@H] flipped to [C@H] for D; extra R-groups preserved as-is.
MANUAL_D_NME = [
    # (chuckles, abbr, name, m_type, m_subtype, rg, ct)

    # NMe-D-Thr: flip alpha only, keep beta stereo
    ("[1*]N(C)[C@H](C([2*])=O)[C@@H](C)O[3*]",
     "D_meT", "NMe-D-threonine", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](CCSC)C([2*])=O",
     "D_meM", "NMe-D-methionine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](CS[3*])C([2*])=O",
     "D_meC", "NMe-D-cysteine", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cncn1[3*])C([2*])=O",
     "D_meH", "NMe-D-histidine", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](CC([3*])=O)C([2*])=O",
     "D_meD", "NMe-D-aspartate", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](CCC([3*])=O)C([2*])=O",
     "D_meE", "NMe-D-glutamate", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](CC(=O)N([3*])[4*])C([2*])=O",
     "D_meN", "NMe-D-asparagine", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CCC(=O)N([3*])[4*])C([2*])=O",
     "D_meQ", "NMe-D-glutamine", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CCCCN([3*])[4*])C([2*])=O",
     "D_meK", "NMe-D-lysine", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CCCN([4*])C(=N)N([3*])[5*])C([2*])=O",
     "D_meR", "NMe-D-arginine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None", "1:backbone_n,2:backbone_c"),

    ("[1*]N(C)[C@H](CC1CCCCC1)C([2*])=O",
     "D_meCha", "NMe-D-cyclohexylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(-c2ccccc2)cc1)C([2*])=O",
     "D_meBip", "NMe-D-biphenylalanine (4)", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1cccc2ccccc12)C([2*])=O",
     "D_me1Nal", "NMe-D-1-naphthylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc2ccccc2c1)C([2*])=O",
     "D_me2Nal", "NMe-D-2-naphthylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # meChg L-form has unusual arm order: [C@H](C([2*])=O)cHex -> flip to [C@@H]
    ("[1*]N(C)[C@@H](C([2*])=O)C1CCCCC1",
     "D_meChg", "NMe-D-cyclohexylglycine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](CC(c1ccccc1)c1ccccc1)C([2*])=O",
     "D_meDip", "NMe-D-diphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](CCC)C([2*])=O",
     "D_meNva", "NMe-D-norvaline", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](CC)C([2*])=O",
     "D_meAbu", "NMe-D-2-aminobutyric acid", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1cccnc1)C([2*])=O",
     "D_me3Pal", "NMe-D-3-pyridylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccncc1)C([2*])=O",
     "D_me4Pal", "NMe-D-4-pyridylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # Phe-variant D-NMe AAs
    ("[1*]N(C)[C@H](Cc1cccc(F)c1)C([2*])=O",
     "D_mePhe_3F", "NMe-D-3-fluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1cccc(OC)c1)C([2*])=O",
     "D_mePhe_3OMe", "NMe-D-3-methoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(Br)cc1)C([2*])=O",
     "D_mePhe_4Br", "NMe-D-4-bromophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(C(F)(F)F)cc1)C([2*])=O",
     "D_mePhe_4CF3", "NMe-D-4-(trifluoromethyl)phenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(Cl)cc1)C([2*])=O",
     "D_mePhe_4Cl", "NMe-D-4-chlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(F)cc1)C([2*])=O",
     "D_mePhe_4F", "NMe-D-4-fluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(OC)cc1)C([2*])=O",
     "D_mePhe_4OMe", "NMe-D-4-methoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(C(C)(C)C)cc1)C([2*])=O",
     "D_mePhe_4tBu", "NMe-D-4-tert-butylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # Trp variants D-NMe: indole N gets [3*] (same as meW/meTrp pattern)
    ("[1*]N(C)[C@H](Cc1cn([3*])c2cc(F)ccc12)C([2*])=O",
     "D_meTrp_5F", "NMe-D-5-fluorotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cn([3*])c2cc(C)ccc12)C([2*])=O",
     "D_meTrp_5Me", "NMe-D-5-methyltryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cn([3*])c2ncccc12)C([2*])=O",
     "D_meTrp_7Aza", "NMe-D-7-azatryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    # NMe-Pro (N-methyl-L-proline): ring N has methyl, no [3*] needed
    ("[1*]N1(C)CCC[C@@H]1C([2*])=O",
     "mePro", "NMe-L-proline", "aa", "modified", _PRO_RG, _PRO_CT),

    ("[1*]N1(C)CCC[C@H]1C([2*])=O",
     "D_mePro", "NMe-D-proline", "aa", "modified", _PRO_RG, _PRO_CT),
]

# ─── SECTION 2: New NMe-AA series ────────────────────────────────────────────
NEW_NME = [
    # NMe-Orn/Dab/Dap: follow meK/meR convention with [3*],[4*] on amine
    ("[1*]N(C)[C@@H](CCCN([3*])[4*])C([2*])=O",
     "meOrn", "NMe-L-ornithine", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CCCN([3*])[4*])C([2*])=O",
     "D_meOrn", "NMe-D-ornithine", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@@H](CCN([3*])[4*])C([2*])=O",
     "meDab", "NMe-L-2,4-diaminobutyric acid", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CCN([3*])[4*])C([2*])=O",
     "D_meDab", "NMe-D-2,4-diaminobutyric acid", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@@H](CN([3*])[4*])C([2*])=O",
     "meDap", "NMe-L-2,3-diaminopropionic acid", "aa", "modified", _NME2_RG, _NME2_CT),

    ("[1*]N(C)[C@H](CN([3*])[4*])C([2*])=O",
     "D_meDap", "NMe-D-2,3-diaminopropionic acid", "aa", "modified", _NME2_RG, _NME2_CT),

    # New NMe-Phe variants (ortho and para)
    ("[1*]N(C)[C@@H](Cc1ccccc1F)C([2*])=O",
     "mePhe_2F", "NMe-L-2-fluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccccc1F)C([2*])=O",
     "D_mePhe_2F", "NMe-D-2-fluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccccc1Cl)C([2*])=O",
     "mePhe_2Cl", "NMe-L-2-chlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccccc1Cl)C([2*])=O",
     "D_mePhe_2Cl", "NMe-D-2-chlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(C)cc1)C([2*])=O",
     "mePhe_4Me", "NMe-L-4-methylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(C)cc1)C([2*])=O",
     "D_mePhe_4Me", "NMe-D-4-methylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(C#N)cc1)C([2*])=O",
     "mePhe_4CN", "NMe-L-4-cyanophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(C#N)cc1)C([2*])=O",
     "D_mePhe_4CN", "NMe-D-4-cyanophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc([N+](=O)[O-])cc1)C([2*])=O",
     "mePhe_4NO2", "NMe-L-4-nitrophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc([N+](=O)[O-])cc1)C([2*])=O",
     "D_mePhe_4NO2", "NMe-D-4-nitrophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(N)cc1)C([2*])=O",
     "mePhe_4NH2", "NMe-L-4-aminophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(N)cc1)C([2*])=O",
     "D_mePhe_4NH2", "NMe-D-4-aminophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(SC)cc1)C([2*])=O",
     "mePhe_4SMe", "NMe-L-4-(methylthio)phenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(SC)cc1)C([2*])=O",
     "D_mePhe_4SMe", "NMe-D-4-(methylthio)phenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(I)cc1)C([2*])=O",
     "mePhe_4I", "NMe-L-4-iodophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(I)cc1)C([2*])=O",
     "D_mePhe_4I", "NMe-D-4-iodophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(Cl)cc1F)C([2*])=O",
     "mePhe_34diF", "NMe-L-3,4-difluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(Cl)cc1F)C([2*])=O",
     "D_mePhe_34diF", "NMe-D-3,4-difluorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # NMe-Trp additional variants
    ("[1*]N(C)[C@@H](Cc1cn([3*])c2c(F)ccc12)C([2*])=O",
     "meTrp_4F", "NMe-L-4-fluorotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cn([3*])c2c(F)ccc12)C([2*])=O",
     "D_meTrp_4F", "NMe-D-4-fluorotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@@H](Cc1cn([3*])c2cc(Br)ccc12)C([2*])=O",
     "meTrp_5Br", "NMe-L-5-bromotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cn([3*])c2cc(Br)ccc12)C([2*])=O",
     "D_meTrp_5Br", "NMe-D-5-bromotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    # NMe-Tyr variants
    ("[1*]N(C)[C@@H](Cc1cc(F)ccc1O[3*])C([2*])=O",
     "meTyr_3F", "NMe-L-3-fluorotyrosine", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cc(F)ccc1O[3*])C([2*])=O",
     "D_meTyr_3F", "NMe-D-3-fluorotyrosine", "aa", "modified", _NME1_RG, _NME1_CT),
]

# ─── SECTION 3: D-forms of batch5 special entries ────────────────────────────
D_BATCH5 = [
    ("[1*]N([3*])[C@H](CCCN([4*])C(=N)N(C)C)C([2*])=O",
     "D_ADMA", "D-asymmetric dimethylarginine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:guanidinium"),

    ("[1*]N([3*])[C@H](CC(=O)OCC=C)C([2*])=O",
     "D_Asp_OAll", "D-aspartate allyl ester", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1nc2ccccc2o1)C([2*])=O",
     "D_BenzoxAla", "D-benzoxazol-2-yl-alanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC#N)C([2*])=O",
     "D_CyanoAla", "D-2-amino-3-cyanopropionic acid", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CSCC=C)C([2*])=O",
     "D_Cys_All", "D-S-allylcysteine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCNC(=O)OCC=C)C([2*])=O",
     "D_Dab_Alloc", "D-Dab(Alloc) - 2,4-diaminobutyric acid N4-allyloxycarbonyl", "aa", "modified",
     _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CNC(=O)OCC=C)C([2*])=O",
     "D_Dap_Alloc", "D-Dap(Alloc) - 2,3-diaminopropionic acid N3-allyloxycarbonyl", "aa", "modified",
     _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCC(=O)OCC=C)C([2*])=O",
     "D_Glu_OAll", "D-glutamate allyl ester", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1nccc2ccccc12)C([2*])=O",
     "D_IsoQuiAla", "D-isoquinolin-1-yl-alanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCNC)C([2*])=O",
     "D_Lys_5Me", "D-N-epsilon-methyllysine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCN(C)C)C([2*])=O",
     "D_Lys_5diMe", "D-N-epsilon-dimethyllysine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCNC(=O)OCC=C)C([2*])=O",
     "D_Lys_Alloc", "D-Lys(Alloc) - lysine N-epsilon-allyloxycarbonyl", "aa", "modified",
     _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCN([4*])C(=O)Oc1ccccc1)C([2*])=O",
     "D_Lys_Cbz", "D-Lys(Cbz) - lysine N-epsilon-benzyloxycarbonyl", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amide_nh"),

    ("[1*]N([3*])[C@H](CCCCNC(c1ccccc1)(c1ccccc1)c1ccc(C)cc1)C([2*])=O",
     "D_Lys_Mtt", "D-Lys(Mtt) - lysine N-epsilon-4-methyltrityl", "aa", "modified",
     _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCN([5*])C(=N)N([4*])[6*])C([2*])=O",
     "D_Nrg", "D-norarginine", "aa", "modified",
     "[H],[OH],[H],[H],[H],[H]",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:guanidinium,6:amine_secondary"),

    ("[1*]N([3*])[C@H](CCCNC(=O)OCC=C)C([2*])=O",
     "D_Orn_Alloc", "D-Orn(Alloc) - ornithine N-delta-allyloxycarbonyl", "aa", "modified",
     _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCN([4*])C(=O)Oc1ccccc1)C([2*])=O",
     "D_Orn_Cbz", "D-Orn(Cbz) - ornithine N-delta-benzyloxycarbonyl", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amide_nh"),

    ("[1*]N([3*])[C@H](Cc1ccc(OC)c(OC)c1)C([2*])=O",
     "D_Phe_34diOMe", "D-3,4-dimethoxyphenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1ccc(C(C)C)cc1)C([2*])=O",
     "D_Phe_4iPr", "D-4-isopropylphenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # D-Pro analogs (flip alpha = last ring C before C=O)
    ("[1*]N1C[C@H](C)C[C@@H]1C([2*])=O",
     "D_Pro_4Me_c", "D-4-methyl-L-proline (cis)", "aa", "modified", _PRO_RG, _PRO_CT),

    ("[1*]N1C[C@@H](C)C[C@@H]1C([2*])=O",
     "D_Pro_4Me_t", "D-4-methyl-L-proline (trans)", "aa", "modified", _PRO_RG, _PRO_CT),

    ("[1*]N1C[C@H](N([3*])[4*])C[C@@H]1C([2*])=O",
     "D_Pro_4NH2_c", "D-4-amino-L-proline (cis)", "aa", "modified", _PRO2_RG, _PRO2_CT),

    ("[1*]N1C[C@@H](N([3*])[4*])C[C@@H]1C([2*])=O",
     "D_Pro_4NH2_t", "D-4-amino-L-proline (trans)", "aa", "modified", _PRO2_RG, _PRO2_CT),

    ("[1*]N([3*])[C@H](Cn1cccc1)C([2*])=O",
     "D_Pyr1Ala", "D-1-pyrazol-1-yl-alanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](COCC=C)C([2*])=O",
     "D_Ser_All", "D-O-allylserine", "aa", "modified", _AA_RG, _AA_CT),

    # D_Thr_All: flip alpha (was [C@H]), keep beta stereo
    ("[1*]N([3*])[C@@H](C([2*])=O)[C@@H](C)OCC=C",
     "D_Thr_All", "D-O-allylthreonine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cn1ccnn1)C([2*])=O",
     "D_TriazAla", "D-1,2,4-triazol-1-yl-alanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1ccc(OCC=C)cc1)C([2*])=O",
     "D_Tyr_All", "D-O-allyltyrosine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCN([5*])C(=O)N([4*])[6*])C([2*])=O",
     "D_hCit", "D-homocitrulline", "aa", "modified",
     "[H],[OH],[H],[H],[H],[H]",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amide_nh,6:amine_secondary"),
]

# ─── SECTION 4: Missing D-alpha-methyl AAs ───────────────────────────────────
D_AME_MISSING = [
    # D-forms of existing L-aMe: flip [C@@] -> [C@] (for those with stereo)
    # or add explicit [C@@] for those originally with no stereo marker
    ("[1*]N([3*])[C@](C)(CCCN([5*])C(=N)N([4*])[6*])C([2*])=O",
     "D_aMeArg", "D-alpha-methylarginine", "aa", "modified",
     "[H],[OH],[H],[H],[H],[H]",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:guanidinium,6:amine_secondary"),

    ("[1*]N([3*])[C@](C)(CS[4*])C([2*])=O",
     "D_aMeCys", "D-alpha-methylcysteine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:thiol"),

    ("[1*]N([3*])[C@](C)(Cc1cncn1[4*])C([2*])=O",
     "D_aMeHis", "D-alpha-methylhistidine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@](C)(CCCCN([4*])[5*])C([2*])=O",
     "D_aMeLys", "D-alpha-methyllysine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@](C)(CCSC)C([2*])=O",
     "D_aMeMet", "D-alpha-methylmethionine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(CCCC)C([2*])=O",
     "D_aMeNle", "D-alpha-methylnorleucine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(CCC)C([2*])=O",
     "D_aMeNva", "D-alpha-methylnorvaline", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N1CCC[C@]1(C)C([2*])=O",
     "D_aMePro", "D-alpha-methylproline", "aa", "modified", _PRO_RG, ""),

    # D_aMeThr: flip alpha [C@](C) -> [C@@](C), keep beta
    ("[1*]N([3*])[C@@](C)(C([2*])=O)[C@@H](C)O",
     "D_aMeThr", "D-alpha-methylthreonine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(O[4*])c(O[5*])c1)C([2*])=O",
     "D_aMeTyr_3OH", "D-alpha-methyl-3-hydroxytyrosine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None", ""),

    # D-forms of existing L-aMePhe variants
    ("[1*]N([3*])[C@](C)(Cc1ccc(Cl)cc1)C([2*])=O",
     "D_aMePhe_4Cl", "D-alpha-methyl-4-chlorophenylalanine", "aa", "modified",
     _AME_RG, "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@](C)(Cc1ccc(F)cc1)C([2*])=O",
     "D_aMePhe_4F", "D-alpha-methyl-4-fluorophenylalanine", "aa", "modified",
     _AME_RG, "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    ("[1*]N([3*])[C@](C)(Cc1ccc(OC)cc1)C([2*])=O",
     "D_aMePhe_4OMe", "D-alpha-methyl-4-methoxyphenylalanine", "aa", "modified",
     _AME_RG, "1:backbone_n,2:backbone_c,3:backbone_n_mod"),
]

# ─── SECTION 5: New L- and D-alpha-methyl AAs ────────────────────────────────
NEW_AME = [
    # New L-aMe: Orn/Dab/Dap
    ("[1*]N([3*])[C@@](C)(CCCN([4*])[5*])C([2*])=O",
     "aMeOrn", "L-alpha-methylornithine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@](C)(CCCN([4*])[5*])C([2*])=O",
     "D_aMeOrn", "D-alpha-methylornithine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@@](C)(CCN([4*])[5*])C([2*])=O",
     "aMeDab", "L-alpha-methyl-2,4-diaminobutyric acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@](C)(CCN([4*])[5*])C([2*])=O",
     "D_aMeDab", "D-alpha-methyl-2,4-diaminobutyric acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@@](C)(CN([4*])[5*])C([2*])=O",
     "aMeDap", "L-alpha-methyl-2,3-diaminopropionic acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@](C)(CN([4*])[5*])C([2*])=O",
     "D_aMeDap", "D-alpha-methyl-2,3-diaminopropionic acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    # New L/D-aMePhe variants
    ("[1*]N([3*])[C@@](C)(Cc1ccc(Br)cc1)C([2*])=O",
     "aMePhe_4Br", "L-alpha-methyl-4-bromophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(Br)cc1)C([2*])=O",
     "D_aMePhe_4Br", "D-alpha-methyl-4-bromophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(C(F)(F)F)cc1)C([2*])=O",
     "aMePhe_4CF3", "L-alpha-methyl-4-(trifluoromethyl)phenylalanine", "aa", "modified",
     _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(C(F)(F)F)cc1)C([2*])=O",
     "D_aMePhe_4CF3", "D-alpha-methyl-4-(trifluoromethyl)phenylalanine", "aa", "modified",
     _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(C#N)cc1)C([2*])=O",
     "aMePhe_4CN", "L-alpha-methyl-4-cyanophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(C#N)cc1)C([2*])=O",
     "D_aMePhe_4CN", "D-alpha-methyl-4-cyanophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(C)cc1)C([2*])=O",
     "aMePhe_4Me", "L-alpha-methyl-4-methylphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(C)cc1)C([2*])=O",
     "D_aMePhe_4Me", "D-alpha-methyl-4-methylphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(C(C)(C)C)cc1)C([2*])=O",
     "aMePhe_4tBu", "L-alpha-methyl-4-tert-butylphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(C(C)(C)C)cc1)C([2*])=O",
     "D_aMePhe_4tBu", "D-alpha-methyl-4-tert-butylphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1cccc(F)c1)C([2*])=O",
     "aMePhe_3F", "L-alpha-methyl-3-fluorophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1cccc(F)c1)C([2*])=O",
     "D_aMePhe_3F", "D-alpha-methyl-3-fluorophenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    # New L/D aMe-naphthyl, cyclohexyl, biphenyl
    ("[1*]N([3*])[C@@](C)(Cc1cccc2ccccc12)C([2*])=O",
     "aMe1Nal", "L-alpha-methyl-1-naphthylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1cccc2ccccc12)C([2*])=O",
     "D_aMe1Nal", "D-alpha-methyl-1-naphthylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc2ccccc2c1)C([2*])=O",
     "aMe2Nal", "L-alpha-methyl-2-naphthylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc2ccccc2c1)C([2*])=O",
     "D_aMe2Nal", "D-alpha-methyl-2-naphthylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(CC1CCCCC1)C([2*])=O",
     "aMeCha", "L-alpha-methylcyclohexylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(CC1CCCCC1)C([2*])=O",
     "D_aMeCha", "D-alpha-methylcyclohexylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(C1CCCCC1)C([2*])=O",
     "aMeChg", "L-alpha-methylcyclohexylglycine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(C1CCCCC1)C([2*])=O",
     "D_aMeChg", "D-alpha-methylcyclohexylglycine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@@](C)(Cc1ccc(-c2ccccc2)cc1)C([2*])=O",
     "aMeBip", "L-alpha-methyl-4-biphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    ("[1*]N([3*])[C@](C)(Cc1ccc(-c2ccccc2)cc1)C([2*])=O",
     "D_aMeBip", "D-alpha-methyl-4-biphenylalanine", "aa", "modified", _AME_RG, _AME_CT),

    # aMeTrp variants (matching aMeTrp's no-stereo convention for L, explicit for D)
    ("[1*]N([3*])C(C)(Cc1cn([4*])c2cc(F)ccc12)C([2*])=O",
     "aMeTrp_5F", "L-alpha-methyl-5-fluorotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@](C)(Cc1cn([4*])c2cc(F)ccc12)C([2*])=O",
     "D_aMeTrp_5F", "D-alpha-methyl-5-fluorotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])C(C)(Cc1cn([4*])c2ncccc12)C([2*])=O",
     "aMeTrp_7Aza", "L-alpha-methyl-7-azatryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@](C)(Cc1cn([4*])c2ncccc12)C([2*])=O",
     "D_aMeTrp_7Aza", "D-alpha-methyl-7-azatryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),
]

# ─── SECTION 6: New L-beta3-homo AAs ─────────────────────────────────────────
# L-form uses [C@H](CC([2*])=O) pattern (matching L-b3hLeu/b3hTyr)
NEW_B3H_L = [
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccccc1",
     "b3hPhe", "L-beta3-homophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)CCSC",
     "b3hMet", "L-beta3-homomethionine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)CCC",
     "b3hNva", "L-beta3-homonorvaline", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)CC1CCCCC1",
     "b3hCha", "L-beta3-homocyclohexylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)C1CCCCC1",
     "b3hChg", "L-beta3-homocyclohexylglycine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1cccc2ccccc12",
     "b3h1Nal", "L-beta3-homo-1-naphthylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc2ccccc2c1",
     "b3h2Nal", "L-beta3-homo-2-naphthylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # Asn and Gln: no R-groups on terminal amide
    ("[1*]N([3*])[C@H](CC([2*])=O)CC(=O)N",
     "b3hAsn", "L-beta3-homoasparagine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)CCC(=O)N",
     "b3hGln", "L-beta3-homoglutamine", "aa", "modified", _AA_RG, _AA_CT),

    # Arg: guanidinium R-groups (following aMeArg/hArg pattern)
    ("[1*]N([3*])[C@H](CC([2*])=O)CCCN([5*])C(=N)N([4*])[6*]",
     "b3hArg", "L-beta3-homoarginine", "aa", "modified",
     "[H],[OH],[H],[H],[H],[H]",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:guanidinium,6:amine_secondary"),

    # His: imidazole [4*] (following b3hHis pattern from hHis)
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1cncn1[4*]",
     "b3hHis", "L-beta3-homohistidine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    # Trp: indole N gets [4*] (following b3hTrp_5F pattern with n([4*]))
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1cn([4*])c2ccccc12",
     "b3hTrp", "L-beta3-homotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),
]

# ─── SECTION 7: New D-beta3-homo AAs ─────────────────────────────────────────
NEW_B3H_D = [
    # D-forms of existing L-b3h (no D-form present)
    ("[1*]N([3*])[C@H](CCCC)CC([2*])=O",
     "D_b3hNle", "D-beta3-homonorleucine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(O)cc1",
     "D_b3hTyr", "D-beta3-homotyrosine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCCCN[4*])CC([2*])=O",
     "D_b3hLys", "D-beta3-homolysine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CC(=O)O[4*]",
     "D_b3hAsp", "D-beta3-homoaspartate", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@H](CCC(=O)O[4*])CC([2*])=O",
     "D_b3hGlu", "D-beta3-homoglutamate", "aa", "modified",
     "[H],[OH],[H],[OH],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(Br)cc1",
     "D_b3hPhe_4Br", "D-beta3-homo-4-bromophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cn([4*])c2cc(F)ccc12",
     "D_b3hTrp_5F", "D-beta3-homo-5-fluorotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),

    # D-forms of new L-b3h (flip [C@H] -> [C@@H])
    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccccc1",
     "D_b3hPhe", "D-beta3-homophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCSC",
     "D_b3hMet", "D-beta3-homomethionine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCC",
     "D_b3hNva", "D-beta3-homonorvaline", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CC1CCCCC1",
     "D_b3hCha", "D-beta3-homocyclohexylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)C1CCCCC1",
     "D_b3hChg", "D-beta3-homocyclohexylglycine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cccc2ccccc12",
     "D_b3h1Nal", "D-beta3-homo-1-naphthylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc2ccccc2c1",
     "D_b3h2Nal", "D-beta3-homo-2-naphthylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CC(=O)N",
     "D_b3hAsn", "D-beta3-homoasparagine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCC(=O)N",
     "D_b3hGln", "D-beta3-homoglutamine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCCN([5*])C(=N)N([4*])[6*]",
     "D_b3hArg", "D-beta3-homoarginine", "aa", "modified",
     "[H],[OH],[H],[H],[H],[H]",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:guanidinium,6:amine_secondary"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cncn1[4*]",
     "D_b3hHis", "D-beta3-homohistidine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cn([4*])c2ccccc12",
     "D_b3hTrp", "D-beta3-homotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod"),
]

# ─── SECTION 8: Extended b3h-Phe halogenated/substituted series ──────────────
B3H_PHE_EXT = [
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc(C#N)cc1",
     "b3hPhe_4CN", "L-beta3-homo-4-cyanophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(C#N)cc1",
     "D_b3hPhe_4CN", "D-beta3-homo-4-cyanophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc(C(F)(F)F)cc1",
     "b3hPhe_4CF3", "L-beta3-homo-4-(trifluoromethyl)phenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(C(F)(F)F)cc1",
     "D_b3hPhe_4CF3", "D-beta3-homo-4-(trifluoromethyl)phenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc(C)cc1",
     "b3hPhe_4Me", "L-beta3-homo-4-methylphenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(C)cc1",
     "D_b3hPhe_4Me", "D-beta3-homo-4-methylphenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc([N+](=O)[O-])cc1",
     "b3hPhe_4NO2", "L-beta3-homo-4-nitrophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc([N+](=O)[O-])cc1",
     "D_b3hPhe_4NO2", "D-beta3-homo-4-nitrophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc(N)cc1",
     "b3hPhe_4NH2", "L-beta3-homo-4-aminophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(N)cc1",
     "D_b3hPhe_4NH2", "D-beta3-homo-4-aminophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccc(I)cc1",
     "b3hPhe_4I", "L-beta3-homo-4-iodophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccc(I)cc1",
     "D_b3hPhe_4I", "D-beta3-homo-4-iodophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # b3h pyridylalanine
    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1cccnc1",
     "b3h3Pal", "L-beta3-homo-3-pyridylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1cccnc1",
     "D_b3h3Pal", "D-beta3-homo-3-pyridylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC([2*])=O)Cc1ccncc1",
     "b3h4Pal", "L-beta3-homo-4-pyridylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CC([2*])=O)Cc1ccncc1",
     "D_b3h4Pal", "D-beta3-homo-4-pyridylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # b3h Orn/Dab/Dap
    ("[1*]N([3*])[C@H](CC([2*])=O)CCCN([4*])[5*]",
     "b3hOrn", "L-beta3-homoornithine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCCN([4*])[5*]",
     "D_b3hOrn", "D-beta3-homoornithine", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@H](CC([2*])=O)CCN([4*])[5*]",
     "b3hDab", "L-beta3-homo-2,4-diaminobutyric acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CCN([4*])[5*]",
     "D_b3hDab", "D-beta3-homo-2,4-diaminobutyric acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@H](CC([2*])=O)CN([4*])[5*]",
     "b3hDap", "L-beta3-homo-2,3-diaminopropionic acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),

    ("[1*]N([3*])[C@@H](CC([2*])=O)CN([4*])[5*]",
     "D_b3hDap", "D-beta3-homo-2,3-diaminopropionic acid", "aa", "modified",
     "[H],[OH],[H],[H],[H],None",
     "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:amine_primary,5:amine_secondary"),
]

# ─── SECTION 9: Misc, click chemistry, homo-AAs, Pro analogs ─────────────────
MISC_MONOMERS = [
    # Selenomethionine (L and D)
    ("[1*]N([3*])[C@@H](CC[Se]C)C([2*])=O",
     "SeM", "L-selenomethionine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CC[Se]C)C([2*])=O",
     "D_SeM", "D-selenomethionine", "aa", "modified", _AA_RG, _AA_CT),

    # Vinyl glycine (2-amino-3-butenoic acid)
    ("[1*]N([3*])[C@@H](C=C)C([2*])=O",
     "VinGly", "L-vinylglycine (2-amino-3-butenoic acid)", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](C=C)C([2*])=O",
     "D_VinGly", "D-vinylglycine", "aa", "modified", _AA_RG, _AA_CT),

    # Homopropargylglycine (2-amino-5-hexynoic acid, click)
    ("[1*]N([3*])[C@@H](CCC#C)C([2*])=O",
     "HPG", "L-homopropargylglycine (2-amino-5-hexynoic acid)", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCC#C)C([2*])=O",
     "D_HPG", "D-homopropargylglycine", "aa", "modified", _AA_RG, _AA_CT),

    # Propargyl cysteine (S-propargyl-Cys, click)
    ("[1*]N([3*])[C@@H](CSCC#C)C([2*])=O",
     "Cys_Prg", "L-S-propargylcysteine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CSCC#C)C([2*])=O",
     "D_Cys_Prg", "D-S-propargylcysteine", "aa", "modified", _AA_RG, _AA_CT),

    # Beta-haloalanines (electrophilic warheads)
    ("[1*]N([3*])[C@@H](CBr)C([2*])=O",
     "Ala_3Br", "L-beta-bromoalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CBr)C([2*])=O",
     "D_Ala_3Br", "D-beta-bromoalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@@H](CCl)C([2*])=O",
     "Ala_3Cl", "L-beta-chloroalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCl)C([2*])=O",
     "D_Ala_3Cl", "D-beta-chloroalanine", "aa", "modified", _AA_RG, _AA_CT),

    # 3-Iodophenylalanine (radio-iodination handle)
    ("[1*]N([3*])[C@@H](Cc1cccc(I)c1)C([2*])=O",
     "Phe_3I", "L-3-iodophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](Cc1cccc(I)c1)C([2*])=O",
     "D_Phe_3I", "D-3-iodophenylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # Pro 4-acetoxy (capped 4-OH-Pro, SPPS compatible)
    ("[1*]N1C[C@@H](OC(=O)C)C[C@H]1C([2*])=O",
     "Pro_4OAc", "L-trans-4-acetoxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO_RG, _PRO_CT),

    ("[1*]N1C[C@@H](OC(=O)C)C[C@@H]1C([2*])=O",
     "D_Pro_4OAc", "D-trans-4-acetoxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO_RG, _PRO_CT),

    # D-Pro 3-hydroxy variants
    ("[1*]N1C[C@H](O[3*])C[C@@H]1C([2*])=O",
     "D_Pro_3OH_c", "D-cis-3-hydroxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO1_RG, _PRO1_CT),

    ("[1*]N1C[C@@H](O[3*])C[C@@H]1C([2*])=O",
     "D_Pro_3OH_t", "D-trans-3-hydroxypyrrolidine-2-carboxylic acid", "aa", "modified",
     _PRO1_RG, _PRO1_CT),

    # D-Pro_4Me3OH and D-Pro_4Me_c/t
    ("[1*]N1CC(C)C(O[3*])[C@@H]1C([2*])=O",
     "D_Pro_4Me3OH", "D-4-methyl-4-hydroxyproline", "aa", "modified", _PRO1_RG, ""),

    # D_Aze (D-azetidine-2-carboxylic acid)
    ("[1*]N1CC[C@H]1C([2*])=O",
     "D_Aze", "D-azetidine-2-carboxylic acid", "aa", "modified", _PRO_RG, _PRO_CT),

    # homo-Tyr L and D
    ("[1*]N([3*])[C@@H](CCc1ccc(O[4*])cc1)C([2*])=O",
     "hTyr", "L-homotyrosine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:hydroxyl"),

    ("[1*]N([3*])[C@H](CCc1ccc(O[4*])cc1)C([2*])=O",
     "D_hTyr", "D-homotyrosine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:hydroxyl"),

    # D-hTrp (hTrp L exists, D missing)
    ("[1*]N([3*])[C@H](CCc1cn([4*])c2ccccc12)C([2*])=O",
     "D_hTrp", "D-homotryptophan", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh"),

    # D-hHis
    ("[1*]N([3*])[C@H](CCc1cn([4*])cn1)C([2*])=O",
     "D_hHis", "D-homohistidine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", ""),

    # homo-Thr L and D (extra CH2 vs Thr)
    ("[1*]N([3*])[C@@H](C[C@@H](C)O[4*])C([2*])=O",
     "hThr", "L-homothreonine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:hydroxyl"),

    ("[1*]N([3*])[C@H](C[C@@H](C)O[4*])C([2*])=O",
     "D_hThr", "D-homothreonine", "aa", "modified",
     "[H],[OH],[H],[H],None,None", "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:hydroxyl"),

    # homo-Cha L and D (extra CH2 vs Cha)
    ("[1*]N([3*])[C@@H](CCC1CCCCC1)C([2*])=O",
     "hCha", "L-homo-cyclohexylalanine", "aa", "modified", _AA_RG, _AA_CT),

    ("[1*]N([3*])[C@H](CCC1CCCCC1)C([2*])=O",
     "D_hCha", "D-homo-cyclohexylalanine", "aa", "modified", _AA_RG, _AA_CT),

    # beta-aminobutyric acid (3-aminobutyric acid, not alpha-AA)
    ("[1*]N([3*])[C@@H](C)CC([2*])=O",
     "bAbu", "L-3-aminobutyric acid (beta-aminobutyric acid)", "aa", "modified", _AA_RG, _AA_CT),
]

# ─── SECTION 10: Additional RCM staple AAs ───────────────────────────────────
# Follows existing S3/S5/S8/R3/R5/R8 pattern
STAPLE_EXTRA = [
    # S4: (S)-alpha-Me with 3-carbon chain to terminal alkene (i,i+3 or i,i+4 depending on ring size)
    ("[1*]N([3*])[C@@](C)(CCC([4*])=C)C([2*])=O",
     "S4", "(S)-alpha-methyl-4-pentenylglycine (RCM staple)", "aa", "modified",
     _STPL_RG, _STPL_CT),

    ("[1*]N([3*])[C@](C)(CCC([4*])=C)C([2*])=O",
     "R4", "(R)-alpha-methyl-4-pentenylglycine (RCM staple)", "aa", "modified",
     _STPL_RG, _STPL_CT),

    # S6: (S)-alpha-Me with 5-carbon chain (i,i+7 ring)
    ("[1*]N([3*])[C@@](C)(CCCCC([4*])=C)C([2*])=O",
     "S6", "(S)-alpha-methyl-6-heptenylglycine (RCM staple)", "aa", "modified",
     _STPL_RG, _STPL_CT),

    ("[1*]N([3*])[C@](C)(CCCCC([4*])=C)C([2*])=O",
     "R6", "(R)-alpha-methyl-6-heptenylglycine (RCM staple)", "aa", "modified",
     _STPL_RG, _STPL_CT),
]

# ─── SECTION 11: Additional acyl cap monomers ────────────────────────────────
CAP_EXTRA = [
    ("[2*]C(=O)CCCCCCC",
     "OctCap", "Octanoyl cap (C8)", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)CCCCCCCCC",
     "DecCap", "Decanoyl cap (C10)", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)CC(C)C",
     "IsoButCap", "3-methylbutanoyl (isovaleryl) cap", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)C(C)(C)C",
     "PivCap", "Pivaloyl (trimethylacetyl) cap", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)C1CCCCC1",
     "CyclohexCap", "Cyclohexanecarbonyl cap", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)Cc1ccccc1",
     "BnCap", "Phenylacetyl (benzyl carbonyl) cap", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)c1ccc2ccccc2c1",
     "NaphtCap", "2-Naphthoyl cap", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)CCCCCCCCCCCCC",
     "MyrCap", "Myristoyl cap (C14)", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)CCCCCCCCCCCCCCC",
     "PalCap", "Palmitoyl cap (C16)", "cap", "cap", _CAP_RG, _CAP_CT),

    ("[2*]C(=O)c1ccccc1",
     "BzCap", "Benzoyl cap", "cap", "cap", _CAP_RG, _CAP_CT),
]

# ─── SECTION 12: More NMe-Phe and extended series ────────────────────────────
MORE_NME_EXT = [
    ("[1*]N(C)[C@@H](Cc1ccc(Cl)cc1F)C([2*])=O",
     "mePhe_3Cl", "NMe-L-3-chlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(Cl)cc1F)C([2*])=O",
     "D_mePhe_3Cl", "NMe-D-3-chlorophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1cccc(Br)c1)C([2*])=O",
     "mePhe_3Br", "NMe-L-3-bromophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1cccc(Br)c1)C([2*])=O",
     "D_mePhe_3Br", "NMe-D-3-bromophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(OC)cc1OC)C([2*])=O",
     "mePhe_34diOMe", "NMe-L-3,4-dimethoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(OC)cc1OC)C([2*])=O",
     "D_mePhe_34diOMe", "NMe-D-3,4-dimethoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(C(C)C)cc1)C([2*])=O",
     "mePhe_4iPr", "NMe-L-4-isopropylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(C(C)C)cc1)C([2*])=O",
     "D_mePhe_4iPr", "NMe-D-4-isopropylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccccc1OC)C([2*])=O",
     "mePhe_2OMe", "NMe-L-2-methoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccccc1OC)C([2*])=O",
     "D_mePhe_2OMe", "NMe-D-2-methoxyphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccccc1C)C([2*])=O",
     "mePhe_2Me", "NMe-L-2-methylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccccc1C)C([2*])=O",
     "D_mePhe_2Me", "NMe-D-2-methylphenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@@H](Cc1ccc(cc1)N(C)C)C([2*])=O",
     "mePhe_4NMe2", "NMe-L-4-dimethylaminophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    ("[1*]N(C)[C@H](Cc1ccc(cc1)N(C)C)C([2*])=O",
     "D_mePhe_4NMe2", "NMe-D-4-dimethylaminophenylalanine", "aa", "modified", _NME_RG, _NME_CT),

    # NMe-Trp 6-F variant
    ("[1*]N(C)[C@@H](Cc1cn([3*])c2ccc(F)cc12)C([2*])=O",
     "meTrp_6F", "NMe-L-6-fluorotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cn([3*])c2ccc(F)cc12)C([2*])=O",
     "D_meTrp_6F", "NMe-D-6-fluorotryptophan", "aa", "modified", _NME1_RG, _NME1_CT),

    # NMe-Tyr 3-I (for radio-iodination)
    ("[1*]N(C)[C@@H](Cc1cc(I)ccc1O[3*])C([2*])=O",
     "meTyr_3I", "NMe-L-3-iodotyrosine", "aa", "modified", _NME1_RG, _NME1_CT),

    ("[1*]N(C)[C@H](Cc1cc(I)ccc1O[3*])C([2*])=O",
     "D_meTyr_3I", "NMe-D-3-iodotyrosine", "aa", "modified", _NME1_RG, _NME1_CT),
]


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor

    rdDepictor.SetPreferCoordGen(True)

    # Read existing abbrs
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    existing = set()
    for m in suppl:
        if m:
            existing.add(m.GetPropsAsDict().get("m_abbr", ""))
    print(f"Existing monomers: {len(existing)}")

    # Collect all entries across sections
    all_entries = (
        MANUAL_D_NME
        + NEW_NME
        + D_BATCH5
        + D_AME_MISSING
        + NEW_AME
        + NEW_B3H_L
        + NEW_B3H_D
        + B3H_PHE_EXT
        + MISC_MONOMERS
        + STAPLE_EXTRA
        + CAP_EXTRA
        + MORE_NME_EXT
    )

    print(f"Entries to process: {len(all_entries)}")

    added = 0
    skipped = 0
    failed = []

    with open(str(SDF_PATH), "a", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for chuckles, abbr, name, mtype, msubtype, rg, ct in all_entries:
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

    print(f"\nDone. Added {added} (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")

    # Verify
    s = list(Chem.SDMolSupplier(str(SDF_PATH), removeHs=False))
    valid = [m for m in s if m]
    print(f"\nFinal valid monomers: {len(valid)}")


if __name__ == "__main__":
    main()
