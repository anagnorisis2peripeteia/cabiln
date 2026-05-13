#!/usr/bin/env python3
"""Batch 12: reduced-amide polyamine residues for cyclam-type macrocycles.

Adds three monomers that let pyPept express polyamine macrocycles (cyclam,
spermine, bis-cyclam-p-xylene linkers, etc.) as reduced-amide
peptidomimetics. Each "residue" is an N-to-CH2 unit with no carbonyl;
inter-residue bonds use the new `reduced_amide` reaction
(`[backbone_n, backbone_c_red] >> [N][C]`).

  redG2     N-CH2-CH2-              (ethano-diamine repeat, 2-C spacer)
  redG3     N-CH2-CH2-CH2-          (propano-diamine repeat, 3-C spacer)
  pXyl     N-CH2-aryl-CH2-N         (para-xylene bis-CH2 linker between two
                                     polyamine macrocycles, as in CP02627)

Cyclam (1,4,8,11-tetraazacyclotetradecane) is then expressible as:
    !1-redG3-redG2-redG3-redG2-!1     (head-to-tail reduced-amide cyclisation)
Bis-cyclam-p-xylene (CP02627) uses two such macrocycles joined by pXyl
through a sidechain (R3) attachment to one N of each cyclam.

Run from repo root:
    python tools/add_monomers_batch12.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH12 = [
    # redG2 — ethano-diamine repeat unit; N-CH2-CH2- (2 C spacer between N's).
    # Backbone: R1=N, R2=terminal CH2 (no carbonyl) tagged backbone_c_red.
    (
        "[1*]N([3*])CC[2*]",
        "redG2",
        "ethano-diamine reduced-amide repeat (-N-CH2-CH2-)",
        "aa", "modified",
        "[H],[H],[H],None,None,None",
        "1:backbone_n,2:backbone_c_red,3:backbone_n_mod",
    ),
    # redG3 — propano-diamine repeat unit; N-CH2-CH2-CH2-.
    (
        "[1*]N([3*])CCC[2*]",
        "redG3",
        "propano-diamine reduced-amide repeat (-N-CH2-CH2-CH2-)",
        "aa", "modified",
        "[H],[H],[H],None,None,None",
        "1:backbone_n,2:backbone_c_red,3:backbone_n_mod",
    ),
    # pXyl — para-xylene bridge: N-CH2-c6h4-CH2-N. Both ends tagged
    # backbone_c_red so each side forms a reduced-amide bond to a cyclam N.
    # Slot 1 carries the "incoming" methylene (treated as the N-terminal
    # endpoint of this linker; the parser writes both sides as reduced amides).
    (
        "[1*]Cc1ccc(C[2*])cc1",
        "pXyl",
        "para-xylene (bis-CH2) reduced-amide linker",
        "aa", "modified",
        "[H],[H],None,None,None,None",
        "1:backbone_c_red,2:backbone_c_red",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH12:
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
