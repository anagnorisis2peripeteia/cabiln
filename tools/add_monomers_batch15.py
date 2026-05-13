#!/usr/bin/env python3
"""Batch 15: phosphine-oxide 3-arm scaffold for CP01557.

Adds PhosOxScaffold — a triaryl-phosphine-oxide scaffold linker (analogous
to TBMB but with P(=O) as the 3-arm hub instead of trisubst-benzene).
Used in CP01557 (aryl-P(=O) hub bridging 3 imidazole-CH2-amide arms).

Structure: each of P(=O)'s 3 aryl substituents is a phenylene that carries
a meta-positioned C(=O) bearing the crosslink R-group. The arm's α-C=O
forms an amide with the αN of the adjacent peptide residue via the new
`aryl_amide_to_backbone_n` reaction.

  PhosOxScaffold   tris(2-carboxyphenyl)-phosphine oxide (3-arm scaffold)
                   slots 4/5/6 expose aryl-amide-C anchors (leaving group [OH]).

Run from repo root:
    python tools/add_monomers_batch15.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH15 = [
    (
        "[4*]C(=O)c1cccc(c1)P(=O)(c1cccc(c1)C([5*])=O)c1cccc(c1)C([6*])=O",
        "PhosOxScaffold",
        "tris(2-carboxyphenyl)-phosphine oxide (3-arm scaffold)",
        "linker", "crosslink",
        "None,None,None,[OH],[OH],[OH]",
        "4:aryl_amide_c,5:aryl_amide_c,6:aryl_amide_c",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH15:
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
