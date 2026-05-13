#!/usr/bin/env python3
"""Batch 14: phosphine-bridged peptidomimetic monomers for CP01557.

Adds two specialised monomers needed to express the cyclic hexapeptidomimetic
CP01557 (an N-methylimidazole-and-aryl-P(=O)-aryl natural-product mimic):

  PhosBridge   Aryl-P(=O)-aryl-amide bridge: an "amide-anchored linker
               residue" — αN at one end, αC=O at the other, with
               -C(=O)-aryl-P(=O)-aryl- as the backbone body. Used as the
               2 inter-residue bridges in CP01557. Bounded composite
               (~11 atoms); 1-of-1 chemistry for this scaffold.

  ImzAla       Imidazole-α residue with a CH2 spacer (β-amino-acid-like):
               -NH-CH2-imidazole(N-Me)-C(=O)-. αN at R1, αC=O at R2,
               imidazole bears the N-methyl observed in CP01557. Reusable
               for other N-methyl-imidazole-containing natural products.

Run from repo root:
    python tools/add_monomers_batch14.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH14 = [
    # PhosBridge — amide-aryl-P(=O)-aryl-amide span. R1=αN, R2=αC=O.
    # Body: -C(=O)-aryl-P(=O)-aryl-C(=O)-  (with the αN at the left end).
    (
        "[1*]NC(=O)c1cccc(c1)P(=O)(c1cccc(c1)C([2*])=O)",
        "PhosBridge",
        "aryl-P(=O)-aryl-amide bridge (CP01557)",
        "aa", "modified",
        "[H],[OH],None,None,None,None",
        "1:backbone_n,2:backbone_c",
    ),
    # ImzAla — N-methyl-imidazole-α residue with CH2 spacer (β-aa-like).
    # Topology: NH(R1)-CH2-imidazole(N-Me)-C(=O)(R2).
    (
        "[1*]N([3*])Cc1nc(C([2*])=O)n(C)c1",
        "ImzAla",
        "N-methyl-imidazol-4-yl-(β)-Gly (CP01557 imidazole-CH2-amide unit)",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH14:
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

    s = list(Chem.SDMolSupplier(str(SDF_PATH), removeHs=False))
    valid = [m for m in s if m]
    print(f"Final valid monomers: {len(valid)}")


if __name__ == "__main__":
    main()
