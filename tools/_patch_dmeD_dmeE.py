#!/usr/bin/env python3
"""Fix m_Rgroups for D_meD and D_meE: slot 3 should be [OH] (carboxyl), not [H]."""
import sys, pathlib, shutil
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

SDF_PATH = pathlib.Path(__file__).parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"
FIX = {
    "D_meD": "[H],[OH],[OH],None,None,None",
    "D_meE": "[H],[OH],[OH],None,None,None",
}

from rdkit import Chem
from rdkit.Chem import SDWriter

suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
all_mols = [(m, m.GetPropsAsDict()) for m in suppl if m is not None]
print(f"Loaded {len(all_mols)}")

tmp = str(SDF_PATH) + ".tmp"
with open(tmp, "w", encoding="utf-8") as fh:
    writer = SDWriter(fh)
    writer.SetKekulize(False)
    for mol, props in all_mols:
        abbr = props.get("m_abbr", "")
        if abbr in FIX:
            mol.SetProp("m_Rgroups", FIX[abbr])
            print(f"  Fixed {abbr}: m_Rgroups -> {FIX[abbr]}")
        writer.write(mol)
    writer.close()

shutil.move(tmp, str(SDF_PATH))
valid = [m for m in Chem.SDMolSupplier(str(SDF_PATH), removeHs=False) if m]
print(f"Done. Final count: {len(valid)}")
