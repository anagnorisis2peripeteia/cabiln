#!/usr/bin/env python3
"""Batch 13: pendant α,β-dehydro-amino-acid residues for chondramide tail.

Adds two unusual α-residues whose α-C is sp2 (carries an exocyclic C=C
double bond). These appear in the pendant DhPro–enamine–bis-acid tail of
chondramide A/B/D macrocycles (CP01144, CP01626, CP02049). Standard
backbone_n / backbone_c chem types — no new reaction needed; only the
structural shape is unusual.

  EnamEt    alpha,beta-dehydro-(2-methyl-but-1-en-1-yl)-glycine
            (α-C=C-CH(Me)-CH2-CH3 side chain) — internal pendant residue.

  EnamDA   alpha,beta-dehydro-(2-carboxyethenyl)-glycine, terminal C-acid
            (α-C=CH-COOH side chain; α-C(=O)OH free) — C-terminal cap.

Run from repo root:
    python tools/add_monomers_batch13.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

BATCH13 = [
    # EnamEt — internal alpha,beta-dehydro residue with ethyl-methyl-vinyl side.
    # Free form: H2N-C(=C(CH3)CH2CH3)-C(=O)-OH. Both R1 and R2 active.
    (
        "[1*]N([3*])/C(=C(\\C)CC)C([2*])=O",
        "EnamEt",
        "alpha,beta-dehydro-(2-methylbut-1-en-1-yl)-glycine",
        "aa", "modified",
        "[H],[OH],[H],None,None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod",
    ),
    # EnamDA — terminal alpha,beta-dehydro residue with vinyl-carboxyl side and free
    # alpha-COOH. C-terminal cap (no R2).
    # Free form: H2N-C(=CH-COOH)-C(=O)OH.
    (
        "[1*]N([3*])/C(=C/C(=O)O)C(=O)O",
        "EnamDA",
        "alpha,beta-dehydro-(2-carboxyethenyl)-glycine (terminal, free acid)",
        "aa", "modified",
        "[H],None,[H],None,None,None",
        "1:backbone_n,3:backbone_n_mod",
    ),
    # _EnamDA — same residue declared as C-terminal cap so the parser picks
    # it up via _s2c_get_cap_libs (which filters on m_type=='cap').
    (
        "[1*]N([3*])/C(=C/C(=O)O)C(=O)O",
        "_EnamDA",
        "alpha,beta-dehydro-bis-carboxyl C-terminal cap (chondramide pendant)",
        "cap", "cap",
        "[H],None,[H],None,None,None",
        "1:backbone_n",
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH13:
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
