#!/usr/bin/env python3
"""Patch v1: fix two categories of data issues in monomers.sdf.

 1. Sidechain carboxyl [4*] on O instead of C
    (hAsp/D_hAsp/hGlu/D_hGlu/b3hAsp/D_b3hAsp/b3hGlu/D_b3hGlu)
    The [4*] dummy must be on the carbonyl C, not on the OH oxygen.

 2. NMe-AAs with sidechain reactive groups whose m_chem_types omits the
    sidechain slot (dummy present in CHUCKLES but not declared in CT).

Run from repo root:
    python tools/patch_monomers_v1.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# ── Category 1: [4*] moved from OH-oxygen to carbonyl-C ──────────────────────
# Old SMILES → New SMILES (only the sidechain carboxyl part changes).
# The correction is: C(=O)O[4*]  →  C([4*])=O   (leaving OH out; dummy on C).
CARBOXYL_FIX = {
    "hAsp":    ("[1*]N([3*])[C@@H](CCC([4*])=O)C([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "D_hAsp":  ("[1*]N([3*])[C@H](CCC([4*])=O)C([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "hGlu":    ("[1*]N([3*])[C@@H](CCCC([4*])=O)C([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "D_hGlu":  ("[1*]N([3*])[C@H](CCCC([4*])=O)C([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    # beta3-homo variants — stereocenter at beta carbon; backbone C is CC([2*])=O
    "b3hAsp":  ("[1*]N([3*])[C@H](CC([2*])=O)CC([4*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "D_b3hAsp":("[1*]N([3*])[C@@H](CC([2*])=O)CC([4*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "b3hGlu":  ("[1*]N([3*])[C@@H](CCC([4*])=O)CC([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
    "D_b3hGlu":("[1*]N([3*])[C@H](CCC([4*])=O)CC([2*])=O",
                "[H],[OH],[H],[OH],None,None",
                "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:carboxyl"),
}

# ── Category 2: NMe-AAs missing sidechain m_chem_types entries ───────────────
# Format: abbr → (new_m_chem_types,)   (m_Rgroups already correct)
# slot 3 type inferred from attachment atom.
_NME_TRP_CT  = "1:backbone_n,2:backbone_c,3:aromatic_nh"
_NME_TYR_CT  = "1:backbone_n,2:backbone_c,3:hydroxyl_phenolic"
_NME_CYS_CT  = "1:backbone_n,2:backbone_c,3:thiol"
_NME_THR_CT  = "1:backbone_n,2:backbone_c,3:hydroxyl"
_NME_HIS_CT  = "1:backbone_n,2:backbone_c,3:aromatic_nh"
_NME_ASP_CT  = "1:backbone_n,2:backbone_c,3:carboxyl"
_NME_GLU_CT  = "1:backbone_n,2:backbone_c,3:carboxyl"
_NME_ASN_CT  = "1:backbone_n,2:backbone_c,3:amide_nh,4:amine_primary"
_NME_GLN_CT  = "1:backbone_n,2:backbone_c,3:amide_nh,4:amine_primary"
_NME_LYS_CT  = "1:backbone_n,2:backbone_c,3:amine_primary,4:amine_secondary"
_NME_ARG_CT  = "1:backbone_n,2:backbone_c,3:guanidinium,4:guanidinium,5:guanidinium"
_NME_ORN_CT  = "1:backbone_n,2:backbone_c,3:amine_primary,4:amine_secondary"
_NME_DAB_CT  = "1:backbone_n,2:backbone_c,3:amine_primary,4:amine_secondary"
_NME_DAP_CT  = "1:backbone_n,2:backbone_c,3:amine_primary,4:amine_secondary"

CT_FIX = {
    # NMe-Trp variants (slot 3 = indole NH)
    "meTrp_5F":    _NME_TRP_CT,
    "meTrp_5Me":   _NME_TRP_CT,
    "meTrp_7Aza":  _NME_TRP_CT,
    "meTrp_5Br":   _NME_TRP_CT,
    "meTrp_6F":    _NME_TRP_CT,
    "D_meW":       _NME_TRP_CT,
    "D_meTrp_5F":  _NME_TRP_CT,
    "D_meTrp_5Me": _NME_TRP_CT,
    "D_meTrp_7Aza":_NME_TRP_CT,
    "D_meTrp_5Br": _NME_TRP_CT,
    "D_meTrp_6F":  _NME_TRP_CT,
    # NMe-Tyr variants (slot 3 = phenolic OH)
    "meTyr_3F":    _NME_TYR_CT,
    "meTyr_3I":    _NME_TYR_CT,
    "D_meTyr_3F":  _NME_TYR_CT,
    "D_meTyr_3I":  _NME_TYR_CT,
    # NMe-D-AAs with single sidechain slot 3
    "D_meC":       _NME_CYS_CT,
    "D_meT":       _NME_THR_CT,
    "D_meH":       _NME_HIS_CT,
    "D_meD":       _NME_ASP_CT,
    "D_meE":       _NME_GLU_CT,
    # NMe-D-AAs with slots 3+4
    "D_meN":       _NME_ASN_CT,
    "D_meQ":       _NME_GLN_CT,
    "D_meK":       _NME_LYS_CT,
    # NMe-D-Arg with slots 3+4+5
    "D_meR":       _NME_ARG_CT,
    # NMe-Orn/Dab/Dap
    "meOrn":       _NME_ORN_CT,
    "D_meOrn":     _NME_ORN_CT,
    "meDab":       _NME_DAB_CT,
    "D_meDab":     _NME_DAB_CT,
    "meDap":       _NME_DAP_CT,
    "D_meDap":     _NME_DAP_CT,
    # beta3-homo-Trp (slot 4 = indole NH; slots 1-3 already in CT)
    "b3hTrp":      "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh",
    "b3hTrp_5F":   "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh",
    "D_b3hTrp_5F": "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh",
    # D_aMeTrp (slot 4 = indole NH; aMe-AAs have slots 1-3 for backbone)
    "D_aMeTrp":    "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aromatic_nh",
}


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor

    rdDepictor.SetPreferCoordGen(True)

    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    all_mols = [(m, m.GetPropsAsDict()) for m in suppl if m is not None]
    print(f"Loaded {len(all_mols)} monomers")

    all_targets = set(CARBOXYL_FIX) | set(CT_FIX)
    patched_c = 0
    patched_ct = 0
    failed = []

    import tempfile, os, shutil
    tmp = str(SDF_PATH) + ".tmp"

    with open(tmp, "w", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for mol, props in all_mols:
            abbr = props.get("m_abbr", "")

            if abbr in CARBOXYL_FIX:
                new_smi, new_rg, new_ct = CARBOXYL_FIX[abbr]
                new_mol = Chem.MolFromSmiles(new_smi)
                if new_mol is None:
                    failed.append((abbr, f"bad SMILES: {new_smi}"))
                    writer.write(mol)  # keep old
                    continue
                rdDepictor.Compute2DCoords(new_mol)
                for k, v in props.items():
                    if k == "m_Rgroups":
                        new_mol.SetProp("m_Rgroups", new_rg)
                    elif k == "m_chem_types":
                        new_mol.SetProp("m_chem_types", new_ct)
                    else:
                        new_mol.SetProp(k, str(v))
                writer.write(new_mol)
                patched_c += 1
                print(f"  PATCHED (carboxyl fix): {abbr}")

            elif abbr in CT_FIX:
                mol.SetProp("m_chem_types", CT_FIX[abbr])
                writer.write(mol)
                patched_ct += 1
                print(f"  PATCHED (CT fix): {abbr}")

            else:
                writer.write(mol)

        writer.close()

    shutil.move(tmp, str(SDF_PATH))

    print(f"\nDone. Patched {patched_c} carboxyl-fix + {patched_ct} CT-fix monomers.")
    if failed:
        print(f"\nFAILED: {failed}")

    # Verify final count
    valid = [m for m in Chem.SDMolSupplier(str(SDF_PATH), removeHs=False) if m]
    print(f"Final monomer count: {len(valid)}")


if __name__ == "__main__":
    main()
