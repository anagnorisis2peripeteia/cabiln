#!/usr/bin/env python3
"""Register degenerate cap monomers with both R1 and R2.

Each degenerate cap has both [1*] and [2*] dummies on the attachment atom,
with leaving group [H] for both. The assembly uses whichever side connects
to the backbone; the unused dummy is removed (implicit H added).
"""
import sys
sys.path.insert(0, "src")

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

SDF_PATH = "src/pyPept/data/monomers.sdf"

# Degenerate caps: name, SMILES with both [1*] and [2*], display name
DEGEN_CAPS = [
    ("Bn",   "[1*]C([2*])c1ccccc1",    "Benzyl (degenerate)"),
    ("Me",   "[1*]C([2*])",             "Methyl (degenerate)"),
    ("Et",   "[1*]C([2*])C",            "Ethyl (degenerate)"),
    ("OBn",  "[1*]O([2*])Cc1ccccc1",    "O-Benzyl (degenerate)"),
    ("OMe",  "[1*]O([2*])C",            "O-Methyl (degenerate)"),
    ("NHBn", "[1*]N([2*])Cc1ccccc1",    "NH-Benzyl (degenerate)"),
]


def register_one(abbr, smiles, name, suppl_mols):
    """Create an RDKit mol for the degenerate cap and return it with props."""
    # Check if already registered
    for m in suppl_mols:
        if m and m.GetPropsAsDict().get("m_abbr") == abbr:
            print(f"  {abbr}: already exists, skipping")
            return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  {abbr}: FAILED to parse SMILES {smiles}")
        return None

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)

    # Set isotopes on dummy atoms: [1*] gets isotope 1, [2*] gets isotope 2
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # The isotope is already set from SMILES parsing ([1*]=1, [2*]=2)
            pass

    mol.SetProp("symbol", abbr)
    mol.SetProp("m_abbr", abbr)
    mol.SetProp("m_name", name)
    mol.SetProp("m_type", "cap")
    mol.SetProp("m_subtype", "cap")
    mol.SetProp("m_Rgroups", "[H],[H],None,None,None,None")
    mol.SetProp("m_degenerate", "true")
    mol.SetProp("_Name", abbr)

    print(f"  {abbr}: OK ({Chem.MolToSmiles(mol)})")
    return mol


def main():
    print("Loading existing monomers...")
    suppl = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
    existing = [m for m in suppl if m is not None]
    print(f"  {len(existing)} monomers loaded")

    print("\nRegistering degenerate caps...")
    new_mols = []
    for abbr, smiles, name in DEGEN_CAPS:
        mol = register_one(abbr, smiles, name, existing)
        if mol:
            new_mols.append(mol)

    if not new_mols:
        print("\nNo new monomers to add.")
        return

    print(f"\nAppending {len(new_mols)} new monomers to SDF...")
    writer = Chem.SDWriter(SDF_PATH)
    writer.SetForceV3000(False)
    for m in existing:
        writer.write(m)
    for m in new_mols:
        writer.write(m)
    writer.close()
    print("Done.")

    # Verify
    suppl2 = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
    count = sum(1 for m in suppl2 if m is not None)
    print(f"  SDF now has {count} monomers")


if __name__ == "__main__":
    main()
