#!/usr/bin/env python3
"""Batch 3: fluorescent dye-lysine conjugate monomers.

All are pre-conjugated Lys building blocks where the epsilon-amine is permanently
modified by the dye. Only backbone R-groups are exposed: R1 (alpha-NH, backbone_n)
and R2 (alpha-C=O, backbone_c). The dye sidechain is inert.

Monomers use manual CHUCKLES format: [1*]N([3*])[C@@H](sidechain)C([2*])=O

Run from repo root:
    python tools/add_monomers_batch3_dyes.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# Standard backbone-only R-groups for dye-Lys monomers
RGROUPS_STD = "[H],[OH],[H],None,None,None"
CHEM_TYPES_STD = "1:backbone_n,2:backbone_c,3:backbone_n_mod"

# (chuckles, abbr, name, m_type, m_subtype)
DYE_MONOMERS = [

    # Dansyl-Lys: Ne-dansyl sulfonamide (5-dimethylaminonaphthalene-1-sulfonyl)
    # ex ~340 nm / em ~510 nm
    ("[1*]N([3*])[C@@H](CCCCNS(=O)(=O)c1cccc2c(N(C)C)cccc12)C([2*])=O",
     "DnLys", "N-epsilon-Dansyl-L-lysine", "aa", "modified"),

    # NBD-Lys: 4-(7-nitrobenzofurazan-4-yl)amino  (benzoxadiazole-amine)
    # ex ~470 nm / em ~540 nm
    ("[1*]N([3*])[C@@H](CCCCNc1cc([N+](=O)[O-])c2nonc2c1)C([2*])=O",
     "NBDLys", "N-epsilon-NBD-L-lysine", "aa", "modified"),

    # AMCA-Lys: 7-amino-4-methylcoumarin-3-acetic acid amide
    # ex ~350 nm / em ~440 nm (blue)
    ("[1*]N([3*])[C@@H](CCCCNC(=O)Cc1c(C)c2ccc(N)cc2oc1=O)C([2*])=O",
     "AMCALys", "N-epsilon-AMCA-L-lysine", "aa", "modified"),

    # Coumarin-343-Lys: julolidine-fused coumarin amide
    # ex ~440 nm / em ~490 nm (cyan)
    ("[1*]N([3*])[C@@H](CCCCNC(=O)c1cc2cc3c4c(c2oc1=O)CCCN4CCC3)C([2*])=O",
     "C343Lys", "N-epsilon-Coumarin343-L-lysine", "aa", "modified"),

    # FITC-Lys: fluorescein-5-isothiocyanate thiourea on epsilon-amine
    # ex ~494 nm / em ~519 nm (green)
    ("[1*]N([3*])[C@@H](CCCCNC(=S)Nc1ccc2c(c1)C(=O)OC21c2ccc(O)cc2Oc2ccc(O)cc21)C([2*])=O",
     "FITCLys", "N-epsilon-FITC-L-lysine (thiourea)", "aa", "modified"),

    # BODIPY-FL-Lys: BODIPY 493/503 C5-pentanoic acid amide
    # ex ~493 nm / em ~503 nm (green)
    ("[1*]N([3*])[C@@H](CCCCNC(=O)CCCCC1=[N+]2C(=Cc3c(C)cc(C)n3[B-]2(F)F)C=C1)C([2*])=O",
     "BDPFLLys", "N-epsilon-BODIPY-FL-C5-L-lysine", "aa", "modified"),

    # 5-TAMRA-Lys: 5-carboxytetramethylrhodamine amide (cationic xanthene)
    # ex ~546 nm / em ~579 nm (orange-red)
    ("[1*]N([3*])[C@@H](CCCCNC(=O)c1ccc(-c2c3ccc(=[N+](C)C)cc-3oc3cc(N(C)C)ccc23)cc1)C([2*])=O",
     "TAMRALys", "N-epsilon-5-TAMRA-L-lysine", "aa", "modified"),

    # Cy3-Lys: non-sulfonated Cy3 trimethine cyanine amide (cationic indolium)
    # ex ~550 nm / em ~570 nm (orange)
    (r"[1*]N([3*])[C@@H](CCCCNC(=O)CCCCCN1/C(=C\C=C\C2=[N+](C)c3ccccc3C2(C)C)c2ccccc21)C([2*])=O",
     "Cy3Lys", "N-epsilon-Cy3-L-lysine (non-sulfonated)", "aa", "modified"),
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

    with open(SDF_PATH, "a", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for chuckles, abbr, name, mtype, msubtype in DYE_MONOMERS:
            if abbr in existing:
                skipped += 1
                print(f"  SKIP {abbr} (exists)")
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
                mol.SetProp("m_Rgroups", RGROUPS_STD)
                mol.SetProp("m_chem_types", CHEM_TYPES_STD)
                writer.write(mol)
                added += 1
                print(f"  ADD  {abbr}: {name}")
            except Exception as exc:
                failed.append((abbr, str(exc)[:100]))

        writer.close()

    print(f"\nDone. Added {added} dye monomers (skipped {skipped} existing).")
    if failed:
        print(f"FAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")


if __name__ == "__main__":
    main()
