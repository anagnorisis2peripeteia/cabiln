"""Stereo validation avoids empty RDKit iteration without trusting stale facts."""

import pytest
from rdkit import Chem

from pyPept.structure import require_supported_stereo


def test_empty_group_vector_does_not_invoke_rdkit_iterator(monkeypatch):
    molecule = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    groups = molecule.GetStereoGroups()
    assert len(groups) == 0

    def unexpected_iteration(_groups):
        pytest.fail("Empty RDKit vector iteration is expensive; check its length")

    monkeypatch.setattr(type(groups), "__iter__", unexpected_iteration)
    require_supported_stereo(molecule)
    require_supported_stereo(None)


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_stereo_validation_observes_new_groups_on_the_same_molecule(group):
    molecule = Chem.RWMol(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    require_supported_stereo(molecule)
    tagged = Chem.MolFromSmiles(f"N[C@@H](C)C(=O)O |{group}:1|")
    kind = tagged.GetStereoGroups()[0].GetGroupType()
    molecule.SetStereoGroups([Chem.CreateStereoGroup(kind, molecule, [1])])
    with pytest.raises(ValueError, match="AND/OR stereochemistry"):
        require_supported_stereo(molecule)
