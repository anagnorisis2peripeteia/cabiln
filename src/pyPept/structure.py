"""Compare molecular graphs without normalizing away chemical differences."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem


def parse_template_smiles(smiles):
    """Read numbered templates without losing slots in older RDKit releases.

    RDKit #8906 clears dummy isotopes when reading CXSMILES, including absolute
    stereo groups. The last Python 3.9 wheel predates that fix. Recover only lost
    numbered dummy atoms from RDKit's plain-SMILES parse, retaining CX metadata.
    Concrete dummies also survive SDF storage; the faulty query atoms do not.
    Ordinary SMILES and unaffected RDKit versions need no second parse.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if (
        molecule is None
        or "|" not in smiles
        or not any(
            atom.GetAtomicNum() == 0 and atom.GetIsotope() == 0
            for atom in molecule.GetAtoms()
        )
    ):
        return molecule
    parameters = Chem.SmilesParserParams()
    parameters.allowCXSMILES = False
    plain = Chem.MolFromSmiles(smiles, parameters)
    if plain is None or plain.GetNumAtoms() != molecule.GetNumAtoms():
        raise ValueError("CXSMILES changed the numbered template's atom indexing")
    molecule = Chem.RWMol(molecule)
    for atom, original in zip(molecule.GetAtoms(), plain.GetAtoms()):
        if (
            atom.GetAtomicNum() == original.GetAtomicNum() == 0
            and not atom.GetIsotope()
            and original.GetIsotope()
        ):
            molecule.ReplaceAtom(atom.GetIdx(), Chem.Atom(original), preserveProps=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule.GetMol()


def require_supported_stereo(molecule) -> None:
    """Do not discard relative or mixture stereo groups during SMILES emission."""
    if molecule is None:
        return
    groups = molecule.GetStereoGroups()
    # Iterating RDKit's empty C++ vector can be much slower than reading its
    # length. Most templates have no groups; no chemistry check is needed then.
    if len(groups) and any(
        group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE
        for group in groups
    ):
        raise ValueError(
            "Enhanced AND/OR stereochemistry is not supported. "
            "Supply an individual stereoisomer with absolute configurations."
        )


@dataclass(frozen=True)
class StructureComparison:
    exact: bool
    compatible: bool
    source_smiles: str
    candidate_smiles: str


def compare_structures(source, candidate) -> StructureComparison:
    """Distinguish equality from filling unspecified input stereochemistry.

    Compatibility requires the same connectivity, charges, isotopes, and every
    stereo assignment in ``source``. Only the candidate may add stereo details.
    No tautomer, protonation, or disconnected-component normalization is used.
    """
    require_supported_stereo(source)
    require_supported_stereo(candidate)
    source_smiles = Chem.MolToSmiles(source)
    candidate_smiles = Chem.MolToSmiles(candidate) if candidate is not None else ""
    if source_smiles == candidate_smiles:
        return StructureComparison(True, True, source_smiles, candidate_smiles)
    compatible = False
    if candidate is not None:
        source_flat, candidate_flat = Chem.Mol(source), Chem.Mol(candidate)
        Chem.RemoveStereochemistry(source_flat)
        Chem.RemoveStereochemistry(candidate_flat)
        compatible = Chem.MolToSmiles(source_flat) == Chem.MolToSmiles(
            candidate_flat
        ) and candidate.HasSubstructMatch(source, useChirality=True)
    return StructureComparison(False, compatible, source_smiles, candidate_smiles)
