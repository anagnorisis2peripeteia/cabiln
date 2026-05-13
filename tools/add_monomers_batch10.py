#!/usr/bin/env python3
"""Batch 10: N-methyl-alpha-thioglycine for polysulfide-bridge DKPs.

Adds N-methyl-alpha-thio-glycine residues whose side chain is an extended
thiol (-S-S-H). A standard disulfide crosslink between two such residues
forms a four-sulfur (tetrasulfide) bridge -- the motif found in
CP00897 (cyclic-N,N'-dimethyl-alpha-thio-Gly DKP with S4 bridge).

The alpha-carbon bears the side-chain sulfur directly (no beta-CH2),
distinguishing this from Cys / meC. Both L and D variants are added
because CP00897 has an L,D meso pair.

  meAtgSS    N-methyl-(S)-alpha-thio-Gly with side chain -S-S-[*]
  D_meAtgSS  N-methyl-(R)-alpha-thio-Gly  "        "       "

Run from repo root:
    python tools/add_monomers_batch10.py
"""

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

# Model after meC: m_Rgroups=[H],[OH],[H],None,None,None
# m_chem_types: 1:backbone_n, 2:backbone_c, 3:thiol
_RG = "[H],[OH],[H],None,None,None"
_CT = "1:backbone_n,2:backbone_c,3:thiol"

BATCH10 = [
    # L-config (matches meC's [C@@H])
    (
        "[1*]N(C)[C@@H](SS[3*])C([2*])=O",
        "meAtgSS",
        "N-methyl-L-alpha-thio-Gly (disulfide-extended)",
        "aa", "modified", _RG, _CT,
    ),
    # D-config
    (
        "[1*]N(C)[C@H](SS[3*])C([2*])=O",
        "D_meAtgSS",
        "N-methyl-D-alpha-thio-Gly (disulfide-extended)",
        "aa", "modified", _RG, _CT,
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

        for smiles, abbr, name, mtype, msubtype, rg, ct in BATCH10:
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
