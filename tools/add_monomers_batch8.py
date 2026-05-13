#!/usr/bin/env python3
"""Batch 8: dehydroproline variants for NRP/chondramide coverage.

Adds three dehydroproline residues that appear in cyclic depsipeptide
natural products (chondramides, jasplakinolide analogues, etc.) but were
absent from the library:

  DhPro34   L-3,4-dehydroproline  (Cβ=Cγ double bond; in CP01144 family)
  D_DhPro34 D-3,4-dehydroproline
  DhPro23   2,3-dehydroproline    (Cα=Cβ; in CP02648; Cα is sp2 so no chirality)

Run from repo root:
    python tools/add_monomers_batch8.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# Pro-type: secondary N in ring, no backbone_n_mod slot
_PRO_RG = "[H],[OH],None,None,None,None"
_PRO_CT = "1:backbone_n,2:backbone_c"

BATCH8 = [
    # ── 3,4-Dehydroproline (L and D) ─────────────────────────────────────────
    # Double bond between Cβ(3) and Cγ(4); Cα is sp3.
    # Ring traversal from N: N-Cδ-Cγ=Cβ-Cα-N(close); exo C=O at Cα.
    # Appears in: chondramide A/B/D (CP01144, CP01626, CP02049), jasplakinolide
    # analogues, telomestatin-type polythiazoles.
    (
        "[1*]N1CC=C[C@@H]1C([2*])=O",
        "DhPro34",
        "L-3,4-dehydroproline (Δ3,4-Pro)",
        "aa", "modified", _PRO_RG, _PRO_CT,
    ),
    (
        "[1*]N1CC=C[C@H]1C([2*])=O",
        "D_DhPro34",
        "D-3,4-dehydroproline (Δ3,4-Pro)",
        "aa", "modified", _PRO_RG, _PRO_CT,
    ),
    # ── 2,3-Dehydroproline ────────────────────────────────────────────────────
    # Double bond between Cα(2) and Cβ(3); Cα is sp2 → no chirality.
    # Ring traversal from N: N-Cδ-Cγ-Cβ=Cα-N(close); exo C=O at sp2 Cα.
    # Appears in: chondramide C (CP02648) and related polythiazole NRPs.
    (
        "[1*]N1CCC=C1C([2*])=O",
        "DhPro23",
        "2,3-dehydroproline (Δ2,3-Pro)",
        "aa", "modified", _PRO_RG, _PRO_CT,
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH8:
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
