#!/usr/bin/env python3
"""Add Tz, TCO, Mal, Aoa, Ald, Dmb monomers to monomers.sdf."""

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor
from rdkit.Chem import AllChem

SDF_PATH = "src/pyPept/data/monomers.sdf"

# Each entry: (smiles, abbr, name, m_type, m_subtype, m_Rgroups, m_chem_types)
# m_Rgroups: comma-separated, 6 positions (R1..R6), None = not present
NEW_MONOMERS = [
    (
        "[2*]C(=O)C([4*])c1nncnn1",
        "Tz",
        "Tetrazine-acetic acid cap",
        "cap", "cap",
        "None,[OH],None,[H],None,None",
        "2:backbone_c,4:tetrazine_c",
    ),
    (
        "[2*]C(=O)CC([4*])C1=CCCCCCC1",
        "TCO",
        "trans-Cyclooctene-propionic acid cap",
        "cap", "cap",
        "None,[OH],None,[H],None,None",
        "2:backbone_c,4:tco_c",
    ),
    (
        "[2*]C(=O)CCN1C(=O)C([4*])=CC1=O",
        "Mal",
        "Maleimidopropionic acid cap",
        "cap", "cap",
        "None,[OH],None,[H],None,None",
        "2:backbone_c,4:maleimide_c",
    ),
    (
        "[1*]N([3*])[C@@H](CON[4*])C([2*])=O",
        "Aoa",
        "L-alpha-aminooxyalanine",
        "aa", "modified",
        "[H],[OH],[H],[H],None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aminooxy",
    ),
    (
        "[1*]N([3*])[C@@H](Cc1ccc([4*]C=O)cc1)C([2*])=O",
        "Ald",
        "L-4-formylphenylalanine",
        "aa", "modified",
        "[H],[OH],[H],[H],None,None",
        "1:backbone_n,2:backbone_c,3:backbone_n_mod,4:aldehyde",
    ),
    (
        "[1*]Cc1ccc(OC)cc1OC",
        "Dmb",
        "2,4-Dimethoxybenzyl N-cap",
        "cap", "protecting",
        "[H],None,None,None,None,None",
        "1:carbon",
    ),
]


def build_mol(smiles: str, abbr: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {abbr}: {smiles}")
    rdDepictor.Compute2DCoords(mol)
    return mol


def main():
    # Load existing monomers to check for duplicates
    existing = set()
    suppl = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
    mols = []
    for m in suppl:
        if m:
            mols.append(m)
            existing.add(m.GetPropsAsDict().get("m_abbr", ""))

    print(f"Existing monomers: {len(mols)}")

    with open(SDF_PATH, "a") as fh:
        writer = SDWriter(fh)
        writer.SetKekulize(False)
        added = 0
        for smiles, abbr, name, mtype, msubtype, rgroups, chem_types in NEW_MONOMERS:
            if abbr in existing:
                print(f"  SKIP {abbr} — already in SDF")
                continue
            mol = build_mol(smiles, abbr)
            mol.SetProp("m_abbr", abbr)
            mol.SetProp("symbol", abbr)
            mol.SetProp("m_name", name)
            mol.SetProp("m_type", mtype)
            mol.SetProp("m_subtype", msubtype)
            mol.SetProp("m_Rgroups", rgroups)
            mol.SetProp("m_chem_types", chem_types)
            writer.write(mol)
            added += 1
            print(f"  ADD  {abbr}: {smiles}")
        writer.close()

    print(f"Done. Added {added} monomers → {SDF_PATH}")


if __name__ == "__main__":
    main()
