#!/usr/bin/env python3
"""Batch 16: ImzAla residue (N-methyl-imidazole-CH2-amide) for CP01557.

Adds ImzAla — an imidazole-α residue with a methylene spacer:
  -NH-CH2-imidazole(N-Me)-C(=O)-
This is one of the 4 "imidazole-amide" residues found in CP01557. R1 is
the αN (amide N at one end), R2 is the αC=O (aryl-bonded carbonyl at the
other end), with the N-methyl-imidazole ring as the residue body.

Run from repo root:
    python tools/add_monomers_batch16.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH16 = [
    (
        "[1*]N([3*])Cc1nc(C([2*])=O)n(C)c1",
        "ImzAla",
        "N-methyl-imidazol-4-yl-(β)-glycine (CP01557 unit)",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH16:
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
