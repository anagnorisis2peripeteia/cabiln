#!/usr/bin/env python3
"""Batch 20: ValAryl — alpha-iPr-glycine with sp3-C-to-aryl chain exit (CP01557).

ValAryl is the "Val without alpha-C=O" residue: alpha-C sp3 with iPr side
chain and alpha-N at R1, where the chain exit is the alpha-C itself
bonding directly to an aryl carbon (typically imidazole-C2) of an
adjacent ImzScaffold. Chem type sp3_c_anchor at R2 pairs with
aryl_c_anchor of the scaffold via aryl_c_c_bond reaction.

  ValAryl    R1: alpha-N (amide N from prev residue)
             R2: sp3 alpha-C (direct C-C bond to aryl anchor)
             R3: N-methyl mod slot on alpha-N (unused in CP01557)
             Body: alpha-C sp3 with iPr branch

Run from repo root:
    python tools/add_monomers_batch20.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH20 = [
    (
        "[1*]N([3*])[C@H](C(C)C)[2*]",
        "ValAryl",
        "alpha-iPr-glycine with sp3-C-to-aryl chain exit (CP01557 Val)",
        "aa", "modified",
        "[H],[H],[H],None,None,None",
        "1:backbone_n,2:sp3_c_anchor,3:backbone_n_mod",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH20:
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
