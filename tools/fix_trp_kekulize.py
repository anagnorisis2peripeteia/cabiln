#!/usr/bin/env python3
"""Fix kekulize-failing Trp-containing entries from batch5.

These entries were written with [nH] on the indole N but require
n([3*]) or n([4*]) to match the existing Trp pattern in the SDF.

Removes the 6 bad entries and re-adds them with correct CHUCKLES.
"""
import sys, pathlib, shutil, tempfile, os
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BAD_ABBRS = {"meTrp_5F", "meTrp_5Me", "meTrp_7Aza", "D_meW", "D_aMeTrp", "b3hTrp_5F"}

# Corrected CHUCKLES entries.
# NMe-Trp: indole N gets [3*] (same as existing meW)  rgroups [H],[OH],[H],None,None,None
# D-aMeTrp: indole N gets [4*] (same as existing aMeTrp), rgroups [H],[OH],[H],[H],None,None
# b3hTrp_5F: indole N gets [4*] following aMeTrp pattern, rgroups [H],[OH],[H],[H],None,None
FIXED = [
    # (chuckles, abbr, name, m_type, m_subtype, rgroups, chem_types)
    ("[1*]N(C)[C@@H](Cc1cn([3*])c2cc(F)ccc12)C([2*])=O",
     "meTrp_5F",  "NMe-L-5-fluorotryptophan",
     "aa","modified","[H],[OH],[H],None,None,None","1:backbone_n,2:backbone_c"),
    ("[1*]N(C)[C@@H](Cc1cn([3*])c2cc(C)ccc12)C([2*])=O",
     "meTrp_5Me", "NMe-L-5-methyltryptophan",
     "aa","modified","[H],[OH],[H],None,None,None","1:backbone_n,2:backbone_c"),
    ("[1*]N(C)[C@@H](Cc1cn([3*])c2ncccc12)C([2*])=O",
     "meTrp_7Aza","NMe-L-7-azatryptophan",
     "aa","modified","[H],[OH],[H],None,None,None","1:backbone_n,2:backbone_c"),
    ("[1*]N(C)[C@H](Cc1cn([3*])c2ccccc12)C([2*])=O",
     "D_meW",     "NMe-D-tryptophan",
     "aa","modified","[H],[OH],[H],None,None,None","1:backbone_n,2:backbone_c"),
    # D_aMeTrp: backbone N gets [1*][3*], alpha-methyl, indole N gets [4*]
    ("[1*]N([3*])[C@](C)(Cc1cn([4*])c2ccccc12)C([2*])=O",
     "D_aMeTrp",  "D-alpha-methyltryptophan (R)",
     "aa","modified","[H],[OH],[H],[H],None,None","1:backbone_n,2:backbone_c,3:backbone_n_mod"),
    # b3hTrp_5F: beta3 backbone, indole N gets [4*]
    ("[1*]N([3*])[C@@H](Cc1cn([4*])c2cc(F)ccc12)CC([2*])=O",
     "b3hTrp_5F", "beta3-Homo-5-fluorotryptophan (L)",
     "aa","modified","[H],[OH],[H],[H],None,None","1:backbone_n,2:backbone_c,3:backbone_n_mod"),
]


def main():
    from rdkit import Chem
    from rdkit.Chem import SDWriter, rdDepictor

    rdDepictor.SetPreferCoordGen(True)

    # Step 1: Read all entries, filtering out bad ones
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False, sanitize=False)
    kept = []
    removed = []
    for m in suppl:
        if m is None:
            continue
        abbr = m.GetPropsAsDict().get("m_abbr", "")
        if abbr in BAD_ABBRS:
            removed.append(abbr)
        else:
            kept.append(m)

    print(f"Removed {len(removed)} bad entries: {removed}")
    print(f"Keeping {len(kept)} entries")

    # Step 2: Write clean SDF to temp file
    tmp = SDF_PATH.with_suffix(".sdf.tmp")
    with open(str(tmp), "w", encoding="utf-8") as fh:
        w = SDWriter(fh)
        w.SetKekulize(False)
        for m in kept:
            try:
                Chem.SanitizeMol(m)
            except Exception:
                try:
                    Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                except Exception:
                    pass
            w.write(m)
        w.close()
        del w

    # Step 3: Atomically replace
    shutil.move(str(tmp), str(SDF_PATH))
    print("SDF rewritten.")

    # Step 4: Append fixed entries
    added = []
    failed = []
    with open(str(SDF_PATH), "a", encoding="utf-8") as fh:
        w = SDWriter(fh)
        w.SetKekulize(False)
        for chuckles, abbr, name, mtype, msubtype, rg, ct in FIXED:
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
                w.write(mol)
                added.append(abbr)
            except Exception as exc:
                failed.append((abbr, str(exc)[:100]))
        w.close()

    print(f"Re-added {len(added)}: {added}")
    if failed:
        print(f"FAILED: {failed}")

    # Step 5: Verify
    s = list(Chem.SDMolSupplier(str(SDF_PATH), removeHs=False))
    valid = [m for m in s if m]
    print(f"\nFinal valid monomers: {len(valid)}")


if __name__ == "__main__":
    main()
