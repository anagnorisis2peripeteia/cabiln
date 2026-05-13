#!/usr/bin/env python3
"""Batch 11: chondramide-family residues for aryl-ether-crosslinked macrocycles.

Adds three monomers required to express chondramide A/B/D (CP01144,
CP01626, CP02049) as a linear peptide chain closed by a single
aryl-ether side-chain crosslink (chondramide C / CP02648 reuses the
same scaffold with an oxazole pendant — covered separately).

  VinylAla   alpha-isopropenyl-glycine (CH2=C(CH3)-CH< backbone)
             — N-terminal residue of the chondramide macrocycle.

  MeBHTyr    N-methyl-beta-hydroxy-tyrosine, phenol-O exposed as crosslink
             stub (R3 tagged 'aryl_phenol_o') — the "ether donor" residue.
             Connects to XlQuat via the new aryl_ether reaction.

  XlQuat     alpha-residue with quaternary-C(Me)(Et) side chain, the
             quaternary C exposed as crosslink stub at R4 (tagged
             'quat_c_anchor') — the "ether acceptor" residue.

The chondramide macrocycle is then:
    VinylAla - MeBHTyr.!1(3,4) - XlQuat.!1(3,4) - DhPro - ...
where !1 is the aryl-ether crosslink between R3 of MeBHTyr and R4 of XlQuat.

Run from repo root:
    python tools/add_monomers_batch11.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH11 = [
    # VinylAla — alpha-isopropenyl-Gly (alpha-(2-methyl-prop-1-en-1-yl)-Gly).
    # Found at the N-terminal of chondramide A/B/D macrocycles.
    (
        "[1*]N([3*])[C@@H](C(=C)C)C([2*])=O",
        "VinylAla",
        "alpha-isopropenyl-glycine",
        "aa", "modified",
        "[H],[OH],[H],None,None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod",
    ),
    # MeBHTyr — N-methyl-beta-hydroxy-Tyr with phenol-O exposed as R3 crosslink.
    # Side chain: beta-CH(OH)-aryl(para-OH, meta-O-crosslink).
    (
        "[1*]N(C)[C@@H]([C@@H](O)c1ccc(O)c(c1)O[3*])C([2*])=O",
        "MeBHTyr",
        "N-methyl-beta-hydroxy-tyrosine (aryl-O crosslink)",
        "aa", "modified",
        "[H],[OH],[H],None,None,None",
        "1:backbone_n,2:backbone_c,3:aryl_phenol_o",
    ),
    # XlQuat — alpha-residue with quaternary-C(Me)(Et) side chain.
    # The quaternary C exposes R4 for an aryl-ether crosslink back to MeBHTyr.
    (
        "[1*]N([3*])[C@@H]([C@@](C)(CC)[4*])C([2*])=O",
        "XlQuat",
        "alpha-(quaternary-C-Me-Et)-glycine — aryl-ether crosslink anchor",
        "aa", "modified",
        "[H],[OH],[H],[H],None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:quat_c_anchor",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH11:
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
