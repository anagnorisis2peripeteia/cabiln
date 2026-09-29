"""Reaction execution preserves atom owners through its selected orientation."""

import pytest
from rdkit import Chem

from pyPept.interfaces.reaction_library import REACTION_INDEX, run_bond_smirks
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence


def owned(smiles, owner):
    molecule = Chem.MolFromSmiles(smiles)
    for atom in molecule.GetAtoms():
        atom.SetIntProp("_residue_idx", owner)
    return molecule


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("isotopes", [(731, 844), (2, 1)])
def test_junction_ownership_uses_actual_reactant_order(reverse, isotopes):
    nitrogen, carbon = isotopes
    amine = owned(f"[{nitrogen}*]N", 17)
    acid = owned(f"[{carbon}*]C(C)=O", 29)
    fragments = (acid, amine) if reverse else (amine, acid)
    sites = (carbon, nitrogen) if reverse else (nitrogen, carbon)
    product = run_bond_smirks(
        *fragments, *sites, REACTION_INDEX[("backbone_n", "backbone_c")], False
    )
    assert Chem.MolToSmiles(product) == "CC(N)=O"
    assert all(
        atom.GetIntProp("_residue_idx") == (17 if atom.GetAtomicNum() == 7 else 29)
        for atom in product.GetAtoms()
    )
    assert Chem.MolToSmiles(amine) == f"[{nitrogen}*]N"
    assert all(atom.GetIntProp("_residue_idx") == 17 for atom in amine.GetAtoms())


def test_later_reaction_step_inherits_from_the_immediate_product():
    # The second step reverses the product's atom order and maps every atom.
    # Original fragment indices therefore cannot recover the final ownership.
    entry = {
        "id": "two-step-lineage-control",
        "slot_a": 1,
        "slot_b": 2,
        "steps": [
            "[1*][N:1].[2*][C:2] >> [N:1][C:2]",
            "[N:1][C:2] >> [C:2][N:1]",
        ],
    }
    product = run_bond_smirks(
        owned("[101*]N", 5), owned("[202*]C", 8), 101, 202, entry, False
    )
    assert Chem.MolToSmiles(product) == "CN"
    assert [atom.GetIntProp("_residue_idx") for atom in product.GetAtoms()] == [8, 5]


def test_reaction_restores_only_owners_lost_at_the_junction(monkeypatch):
    amine = owned("[101*]NCCCCCC", 17)
    acid = owned("[202*]C(C)=O", 29)
    assigned = []
    set_property = Chem.Atom.SetIntProp

    def record_assignment(atom, key, value):
        if key == "_residue_idx":
            assigned.append((atom.GetSymbol(), value))
        return set_property(atom, key, value)

    monkeypatch.setattr(Chem.Atom, "SetIntProp", record_assignment)
    product = run_bond_smirks(
        amine, acid, 101, 202, REACTION_INDEX[("backbone_n", "backbone_c")], False
    )
    assert Chem.MolToSmiles(product) == "CCCCCCNC(C)=O"
    assert sorted(assigned) == [("C", 29), ("N", 17), ("O", 29)]
    assert sorted(atom.GetIntProp("_residue_idx") for atom in product.GetAtoms()) == (
        [17] * 7 + [29] * 3
    )


@pytest.mark.parametrize(
    "notation",
    [
        "ac-C.!1(4,4)-A-C.!1-am",
        "C.!1(4,4)%C.!1",
        "K.[G(4,2).ac(1,2)]-A",
        "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am",
        "ac-S5.!1(4,4)-A-A-A-S5.!1(4,4)-am",
        "ac-D.!1(4,3)-A.!1(3,4)-G-am",
    ],
)
def test_surviving_atom_tracers_keep_their_occurrence_owner(notation):
    sequence = Sequence(notation)
    for owner, monomer in enumerate(sequence.s_monomers):
        # Definitions can be shared by repeated occurrences within a Sequence.
        monomer["m_romol"] = Chem.Mol(monomer["m_romol"])
        for atom in monomer["m_romol"].GetAtoms():
            if atom.GetAtomicNum():
                atom.SetIsotope(100 + owner)
    product = Molecule(sequence, depiction=None).mol
    assert all(atom.HasProp("_residue_idx") for atom in product.GetAtoms())
    traced = [atom for atom in product.GetAtoms() if atom.GetIsotope() >= 100]
    assert traced
    assert all(
        atom.GetIntProp("_residue_idx") == atom.GetIsotope() - 100 for atom in traced
    )


@pytest.mark.parametrize(
    "notation",
    ["ac-C.!1(4,4)-A-C.!1-am", "C.!1(4,4)%C.!1", "A-G", "K.[G(4,2).ac(1,2)]-A"],
)
def test_assembled_ownership_is_complete_without_neighbor_inference(notation):
    assembly = Molecule(Sequence(notation), depiction=None)
    mapping = assembly.get_residue_atom_map()
    assert sorted(atom for atoms in mapping.values() for atom in atoms) == list(
        range(assembly.mol.GetNumAtoms())
    )
    for owner, indices in mapping.items():
        assert all(
            assembly.mol.GetAtomWithIdx(index).GetIntProp("_residue_idx") == owner
            for index in indices
        )
