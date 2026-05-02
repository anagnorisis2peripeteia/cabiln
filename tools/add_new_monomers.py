#!/usr/bin/env python3
"""Add Tz, TCO, Mal, Aoa, Ald, Dmb monomers to monomers.sdf.

CHUCKLES validation
-------------------
After building each mol, _validate_chuckles() calls infer_chem_type on
every [n*] dummy and asserts the detected type matches the declared
m_chem_types string.  This catches placement errors like the original
Ald bug ([4*]C=O vs C([4*])=O) before they reach the SDF.
"""

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor

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
        "[1*]N([3*])[C@@H](Cc1ccc(C([4*])=O)cc1)C([2*])=O",
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

# Slots that are identified purely by slot number (1-3 backbone); infer_chem_type
# uses position-based logic for these, so we skip SMARTS validation for them.
_BACKBONE_SLOTS = {1, 2, 3}


def _validate_chuckles(mol: Chem.Mol, abbr: str, chem_types_str: str) -> None:
    """Assert every [n*] dummy matches its declared chem_type via infer_chem_type.

    Catches CHUCKLES authoring mistakes like placing [4*] on the wrong atom
    before the mol is written to the SDF.
    """
    from pyPept.interfaces.reaction_library import infer_chem_type

    declared = {}
    for tok in chem_types_str.split(","):
        slot_s, ct = tok.strip().split(":")
        declared[int(slot_s)] = ct.strip()

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        slot = atom.GetIsotope()
        if slot not in declared or slot in _BACKBONE_SLOTS:
            continue
        expected = declared[slot]
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            continue
        attach_idx = neighbors[0].GetIdx()
        detected = infer_chem_type(mol, attach_idx, slot=slot)
        if detected != expected:
            raise ValueError(
                f"{abbr} R{slot}: CHUCKLES validation failed — "
                f"declared '{expected}' but infer_chem_type detected '{detected}'.\n"
                f"  Check dummy placement in SMILES: {Chem.MolToSmiles(mol)}"
            )
    print(f"    validated: {chem_types_str}")


def build_mol(smiles: str, abbr: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {abbr}: {smiles}")
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    return mol


def main():
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
                print(f"  SKIP {abbr} (already in SDF)")
                continue
            mol = build_mol(smiles, abbr)
            _validate_chuckles(mol, abbr, chem_types)   # raises if dummy misplaced
            mol.SetProp("m_abbr",      abbr)
            mol.SetProp("symbol",      abbr)
            mol.SetProp("m_name",      name)
            mol.SetProp("m_type",      mtype)
            mol.SetProp("m_subtype",   msubtype)
            mol.SetProp("m_Rgroups",   rgroups)
            mol.SetProp("m_chem_types", chem_types)
            writer.write(mol)
            added += 1
            print(f"  ADD  {abbr}: {smiles}")
        writer.close()

    print(f"Done. Added {added} monomers.")


if __name__ == "__main__":
    main()
