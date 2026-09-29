"""A compatibility baseline must never certify an altered or custom definition."""

import pytest
from rdkit import Chem

from pyPept.library_quality import (
    audit_monomer,
    definition_hash,
    quality_for_monomer,
    quality_manifest,
)
from pyPept.molecule import Molecule
from pyPept.monomer_store import _load_sdf
from pyPept.sequence import Sequence


def definition(symbol):
    return Chem.Mol(_load_sdf()[1][symbol])


def product(source):
    return Chem.MolToSmiles(Molecule(Sequence(source), depiction=None).mol)


def test_quality_lookup_is_bound_to_definition_and_returns_detached_facts():
    original = definition("meC")
    actual = quality_for_monomer(original)
    assert actual["status"] == "review_required"
    assert actual["activation_status"] == "legacy_slot_match"
    assert "legacy_r3_sidechain" in {issue["code"] for issue in actual["issues"]}
    actual["issues"].clear()
    assert quality_for_monomer(original)["issues"]
    changed = Chem.Mol(original)
    changed.GetAtomWithIdx(0).SetIsotope(13)
    assert quality_for_monomer(changed)["status"] == "unreviewed"
    assert definition_hash(changed) != definition_hash(original)


@pytest.mark.parametrize(
    "field,value",
    [
        ("symbol", "RegisteredMeC"),
        ("m_abbr", "DifferentAlias"),
        ("m_Rgroups", "[H],[OH],[Cl]"),
        ("m_chem_types", "1:backbone_n,2:backbone_c,3:alcohol"),
    ],
)
def test_custom_names_and_changed_attachment_metadata_are_unreviewed(field, value):
    molecule = definition("meC")
    molecule.SetProp(field, value)
    result = quality_for_monomer(molecule)
    assert result["status"] == "unreviewed"
    assert result["issues"] == []
    assert result["activation_status"] is None


def test_lookup_does_not_reactivate_definitions(monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("The palette must only read existing audit facts")

    monkeypatch.setattr("pyPept.interfaces.monomer_pipeline.pre_activate", unexpected)
    result = quality_for_monomer(definition("A"))
    assert result["status"] == "no_known_exception"
    assert result["activation_status"] == "match"


def test_definition_fingerprint_ignores_atom_traversal_and_coordinates():
    molecule = definition("meC")
    reordered = Chem.RenumberAtoms(
        molecule, list(reversed(range(molecule.GetNumAtoms())))
    )
    # RenumberAtoms returns structure only; retain the same named definition.
    for name in molecule.GetPropNames():
        reordered.SetProp(name, molecule.GetProp(name))
    reordered.RemoveAllConformers()
    assert definition_hash(reordered) == definition_hash(molecule)
    assert quality_for_monomer(reordered) == quality_for_monomer(molecule)


def test_manifest_is_detached_and_has_an_explicit_baseline_for_every_record():
    manifest = quality_manifest()
    assert manifest["summary"]["records"] == len(manifest["records"]) == 1128
    assert all(record["definition_hash"] for record in manifest["records"].values())
    manifest["records"].clear()
    assert quality_manifest()["records"]


def test_legacy_r3_has_its_original_thiol_product_and_connection():
    expected = Chem.MolFromSmiles("CN[C@@H](CS)C(=O)O")
    assert product("meC") == Chem.MolToSmiles(expected)
    disulfide = Chem.MolFromSmiles("CN[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O")
    assert product("meC.[C(3,4)]") == Chem.MolToSmiles(disulfide)
    record = audit_monomer(definition("meC"))
    assert record["activation"]["status"] == "legacy_slot_match"
    slots = {
        atom.GetIsotope()
        for atom in definition("meC").GetAtoms()
        if atom.GetAtomicNum() == 0
    }
    assert slots == {1, 2, 3}


def test_sparse_terminal_metadata_is_visible_without_changing_restoration():
    result = quality_for_monomer(definition("Mpa"))
    sparse = next(
        issue
        for issue in result["issues"]
        if issue["code"] == "missing_leaving_metadata"
    )
    assert sparse["slots"] == [3]
    expected = Chem.MolFromSmiles("OC(=O)[C@H](C)S")
    assert product("Mpa") == Chem.MolToSmiles(expected)
    disulfide = Chem.MolFromSmiles("OC(=O)[C@H](C)SSC[C@H](N)C(=O)O")
    assert product("Mpa.[C(3,4)]") == Chem.MolToSmiles(disulfide)


def test_unspecified_stereo_is_disclosed_without_inventing_a_configuration():
    molecule = definition("aMeLeu")
    result = quality_for_monomer(molecule)
    assert any(
        issue["code"] == "unspecified_tetrahedral_stereo" for issue in result["issues"]
    )
    expected = Chem.MolFromSmiles("CC(C)CC(C)(N)C(=O)O")
    assert product("aMeLeu") == Chem.MolToSmiles(expected)
    assigned = Chem.MolFromSmiles("CC(C)C[C@](C)(N)C(=O)O")
    assert product("aMeLeu") != Chem.MolToSmiles(assigned)


@pytest.mark.parametrize(
    "symbol,expected_smiles",
    [
        ("Hsl", "N[C@H]1CCOC1=O"),
        ("TATA", "CCC(=O)N1CN(C(=O)CC)CN(C(=O)CC)C1"),
    ],
)
def test_activation_failures_keep_the_existing_standalone_product(
    symbol, expected_smiles
):
    result = quality_for_monomer(definition(symbol))
    assert result["activation_status"] == "activation_failed"
    assert any(issue["code"] == "activation_exception" for issue in result["issues"])
    assert product(symbol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected_smiles))


def test_declared_chemistry_difference_does_not_change_lysine_acylation():
    result = quality_for_monomer(definition("K"))
    issue = next(
        item
        for item in result["issues"]
        if item["code"] == "chemistry_declaration_difference"
    )
    assert issue["sites"] == [
        {"slot": 5, "declared": "amine_secondary", "effective": "amine_primary"}
    ]
    expected = Chem.MolFromSmiles("N[C@@H](CCCCNC(C)=O)C(=O)O")
    assert product("K.Ac(4,2)") == Chem.MolToSmiles(expected)
