#!/usr/bin/env python3
"""Batch 17: ImzIprVal — alpha-iPr-imidazol-2-yl-glycine for CP01557.

Adds ImzIprVal, the chemistry-honest single residue covering one full
"Val + imidazole-via-C2 + αC=O on imidazole-C4" unit of CP01557. The Val
alpha-C is bonded directly to the imidazole C2 (no methylene spacer —
unlike histidine), and the imidazole is N1-methylated and carries the
forward αC=O at its C4.

This is one residue per imidazole-junction (NOT a composite of separable
residues): the imidazole-C2-Val sp3-aryl C-C bond is intrinsic to the
biosynthetic unit, similar to how dehydrohistidine has its α-C bonded to
the imidazole side chain. The C5-CH2-amide branch on the same imidazole
will be handled by adjacent residues in the CP01557 chain.

  ImzIprVal    R1: αN (amide-N of Val)
               R2: αC=O on imidazole-C4 (forward amide to next residue)
               R3: N-methyl mod slot on αN (unused in CP01557)
               Body: αC sp3 with iPr + N-methyl-imidazole

Run from repo root:
    python tools/add_monomers_batch17.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH17 = [
    (
        "[1*]N([3*])[C@H](C(C)C)c1nc(C([2*])=O)n(C)c1",
        "ImzIprVal",
        "alpha-iPr-(1-methylimidazol-2-yl)-glycine (CP01557 unit)",
        "aa", "modified",
        "[H],[OH],[H],None,None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod",
    ),
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH17:
            if abbr in existing:
                skipped += 1
                print(f"  SKIP {abbr} (already present)")
                continue
            try:
                mol = Chem.MolFromSmiles(smiles)
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
                print(f"  ADD  {abbr}")
            except Exception as exc:
                failed.append((abbr, str(exc)[:120]))
                print(f"  FAIL {abbr}: {exc}")

        writer.close()

    print(f"\nDone. Added {added} (skipped {skipped} existing).")
    if failed:
        print(f"\nFAILED ({len(failed)}):")
        for abbr, err in failed:
            print(f"  {abbr}: {err}")


if __name__ == "__main__":
    main()
