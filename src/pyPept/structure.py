"""Compare molecular graphs without normalizing away chemical differences."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem


def require_supported_stereo(molecule) -> None:
    """Do not discard relative or mixture stereo groups during SMILES emission."""
    if molecule is not None and any(
        group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE
        for group in molecule.GetStereoGroups()
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
