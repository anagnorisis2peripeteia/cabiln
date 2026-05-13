#!/usr/bin/env python3
"""Batch 9: α-hydroxy acids for cyclodepsipeptide coverage.

Adds the two most common α-hydroxy acid building blocks found in natural
product cyclodepsipeptides, each as a D/L pair to match the existing
D_Lac/L_Lac convention. All four use the ester backbone
(`1:backbone_o,2:backbone_c`) introduced for lactic / glycolic acid.

  Hiv / D_Hiv     α-hydroxyisovaleric acid (2-hydroxy-3-methylbutanoic)
                  — beauvericin, enniatins A/B, bassianolide, valinomycin
  PhLac / D_PhLac 3-phenyllactic acid (2-hydroxy-3-phenylpropanoic)
                  — PF1022, didemnin, AM-toxin, bassianolide variants

Run from repo root:
    python tools/add_monomers_batch9.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# Ester backbone: O-attached N-side, C(=O)-attached C-side.
_DEPSI_RG = "[H],[OH]"
_DEPSI_CT = "1:backbone_o,2:backbone_c"

BATCH9 = [
    # ── α-Hydroxyisovaleric acid (Hiv) ───────────────────────────────────────
    # (S)-2-hydroxy-3-methylbutanoic acid = L-Hiv (matches L_Lac chirality)
    (
        "[1*]O[C@@H](C(C)C)C([2*])=O",
        "Hiv",
        "L-alpha-hydroxyisovaleric acid",
        "aa", "depsi", _DEPSI_RG, _DEPSI_CT,
    ),
    # (R)-2-hydroxy-3-methylbutanoic acid = D-Hiv (beauvericin/enniatin)
    (
        "[1*]O[C@H](C(C)C)C([2*])=O",
        "D_Hiv",
        "D-alpha-hydroxyisovaleric acid",
        "aa", "depsi", _DEPSI_RG, _DEPSI_CT,
    ),
    # ── 3-Phenyllactic acid (PhLac) ──────────────────────────────────────────
    # (S)-2-hydroxy-3-phenylpropanoic acid = L-PhLac
    (
        "[1*]O[C@@H](Cc1ccccc1)C([2*])=O",
        "PhLac",
        "L-3-phenyllactic acid",
        "aa", "depsi", _DEPSI_RG, _DEPSI_CT,
    ),
    # (R)-2-hydroxy-3-phenylpropanoic acid = D-PhLac (PF1022, didemnin)
    (
        "[1*]O[C@H](Cc1ccccc1)C([2*])=O",
        "D_PhLac",
        "D-3-phenyllactic acid",
        "aa", "depsi", _DEPSI_RG, _DEPSI_CT,
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH9:
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
