#!/usr/bin/env python3
"""Register TBMB (1,3,5-Tris(bromomethyl)benzene) crosslinker monomer."""

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor

SDF_PATH = "src/pyPept/data/monomers.sdf"

suppl = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
existing = set()
for m in suppl:
    if m:
        existing.add(m.GetPropsAsDict().get("m_abbr", ""))
print(f"Existing monomers: {len(existing)}")

if "TBMB" in existing:
    print("TBMB already registered — skipping")
    raise SystemExit(0)

smiles = "[4*]C(Br)c1cc(C([5*])Br)cc(C([6*])Br)c1"
mol = Chem.MolFromSmiles(smiles)
assert mol is not None, "Invalid SMILES"
rdDepictor.SetPreferCoordGen(True)
rdDepictor.Compute2DCoords(mol)

from pyPept.interfaces.reaction_library import infer_chem_type

chem_types = "4:alkyl_halide_c,5:alkyl_halide_c,6:alkyl_halide_c"
for atom in mol.GetAtoms():
    if atom.GetAtomicNum() == 0:
        slot = atom.GetIsotope()
        nb = atom.GetNeighbors()[0]
        ct = infer_chem_type(mol, nb.GetIdx(), slot=slot)
        assert ct == "alkyl_halide_c", f"R{slot}: expected alkyl_halide_c, got {ct}"
print("Validation passed: all 3 arms are alkyl_halide_c")

mol.SetProp("m_abbr",       "TBMB")
mol.SetProp("symbol",       "TBMB")
mol.SetProp("m_name",       "1,3,5-Tris(bromomethyl)benzene")
mol.SetProp("m_type",       "linker")
mol.SetProp("m_subtype",    "crosslink")
mol.SetProp("m_Rgroups",    "None,None,None,[H],[H],[H]")
mol.SetProp("m_chem_types", chem_types)

with open(SDF_PATH, "a") as fh:
    writer = SDWriter(fh)
    writer.SetKekulize(False)
    writer.write(mol)
    writer.close()

print("TBMB added to monomers.sdf")
