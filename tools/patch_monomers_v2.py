#!/usr/bin/env python3
"""Patch monomers.sdf: revert sidechain carboxyl from C-attachment to O-attachment.

The carboxyl chem_type now places the dummy on the OH oxygen (O-attachment):
  CHUCKLES: C(=O)O[4*]  (dummy replaces the H on the hydroxyl)
  leaving group: [H]

patch_monomers_v1 previously moved these to C-attachment: C([4*])=O / LG=[OH].
This patch reverts them.  m_chem_types stays unchanged (4:carboxyl still correct).

Run from repo root:
    python tools/patch_monomers_v2.py
"""

import sys, pathlib, shutil
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"

rdDepictor.SetPreferCoordGen(True)

# Format: abbr -> (new_chuckles_smiles, new_m_Rgroups)
# m_chem_types stays — "4:carboxyl" is correct for both attachment forms.
_REVERT = {
    # α-amino acid homo-variants (hAsp/hGlu = homoAsp/homoGlu)
    "hAsp":    ("[1*]N([3*])[C@@H](CCC(=O)O[4*])C([2*])=O",
                "[H],[OH],[H],[H],None,None"),
    "D_hAsp":  ("[1*]N([3*])[C@H](CCC(=O)O[4*])C([2*])=O",
                "[H],[OH],[H],[H],None,None"),
    "hGlu":    ("[1*]N([3*])[C@@H](CCCC(=O)O[4*])C([2*])=O",
                "[H],[OH],[H],[H],None,None"),
    "D_hGlu":  ("[1*]N([3*])[C@H](CCCC(=O)O[4*])C([2*])=O",
                "[H],[OH],[H],[H],None,None"),
    # beta-3-homo variants (backbone C is CC([2*])=O, sidechain CC(=O)O[4*])
    "b3hAsp":  ("[1*]N([3*])[C@H](CC([2*])=O)CC(=O)O[4*]",
                "[H],[OH],[H],[H],None,None"),
    "D_b3hAsp":("[1*]N([3*])[C@@H](CC([2*])=O)CC(=O)O[4*]",
                "[H],[OH],[H],[H],None,None"),
    "b3hGlu":  ("[1*]N([3*])[C@@H](CCC(=O)O[4*])CC([2*])=O",
                "[H],[OH],[H],[H],None,None"),
    "D_b3hGlu":("[1*]N([3*])[C@H](CCC(=O)O[4*])CC([2*])=O",
                "[H],[OH],[H],[H],None,None"),
}


def main():
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
    all_mols = [(m, m.GetPropsAsDict()) for m in suppl if m is not None]
    print(f"Loaded {len(all_mols)} monomers")

    patched = []
    failed = []
    tmp = str(SDF_PATH) + ".tmp"

    with open(tmp, "w", encoding="utf-8") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)

        for mol, props in all_mols:
            abbr = props.get("m_abbr", "")
            if abbr in _REVERT:
                new_smi, new_rg = _REVERT[abbr]
                new_mol = Chem.MolFromSmiles(new_smi)
                if new_mol is None:
                    print(f"  ERROR: bad SMILES for {abbr}: {new_smi}")
                    failed.append(abbr)
                    writer.write(mol)
                    continue
                rdDepictor.Compute2DCoords(new_mol)
                # Copy all properties; override m_Rgroups
                for k, v in props.items():
                    new_mol.SetProp(k, new_rg if k == "m_Rgroups" else str(v))
                writer.write(new_mol)
                patched.append(abbr)
                print(f"  Patched {abbr}: C(=O)O[4*] O-attachment, LG=[H]")
            else:
                writer.write(mol)

        writer.close()

    shutil.move(tmp, str(SDF_PATH))

    valid = [m for m in Chem.SDMolSupplier(str(SDF_PATH), removeHs=False) if m]
    print(f"\nPatched {len(patched)}: {patched}")
    if failed:
        print(f"Failed: {failed}")
    print(f"Final SDF count: {len(valid)}")


if __name__ == "__main__":
    main()
