#!/usr/bin/env python3
"""Batch 18: ImzScaffold — 3-arm N-methyl-imidazole scaffold for CP01557.

CP01557 has 4 imidazole 3-way junctions, each with:
  - C2: bonded to a Val αC (sp3-aryl C-C)
  - C4: bonded to αC=O of a backwards-amide
  - C5: bonded to CH2 of a forwards-amide
  - N1: bears N-methyl
  - N3: free

This is structurally analogous to TBMB (3-arm benzene scaffold). Modelled
as a `linker/crosslink` monomer with 3 anchor slots (R4, R5, R6) on the
three ring carbons. Each anchor is `aryl_c_anchor` (sp2 aryl-C bonding to
another atom). The bonded partners come from adjacent residues:
  - R4 (C2) bonds to Val αC via aryl_c_aryl_c reaction (sp3 partner)
  - R5 (C4) bonds to amide αC=O via aryl_c_amide_c reaction (sp2 partner)
  - R6 (C5) bonds to amide CH2 via aryl_c_methylene_c reaction (sp3 partner)

Run from repo root:
    python tools/add_monomers_batch18.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH18 = [
    (
        "[4*]c1nc([5*])n(C)c1[6*]",
        "ImzScaffold",
        "N-methyl-imidazole 3-arm scaffold (C2/C4/C5 anchors)",
        "linker", "crosslink",
        "None,None,None,[H],[H],[H]",
        "4:aryl_c_anchor,5:aryl_c_anchor,6:aryl_c_anchor",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH18:
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
