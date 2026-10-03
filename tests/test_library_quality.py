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
    original = definition("ImzScaffold")
    actual = quality_for_monomer(original)
    assert actual["status"] == "no_known_exception"
    assert actual["activation_status"] == "match"
    assert "authored_attachment_sites" in {issue["code"] for issue in actual["issues"]}
    assert all(issue["severity"] == "info" for issue in actual["issues"])
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
        ("m_name", "Different named identity"),
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


def test_migrated_thiol_has_its_original_product_and_connection():
    expected = Chem.MolFromSmiles("CN[C@@H](CS)C(=O)O")
    assert product("meC") == Chem.MolToSmiles(expected)
    disulfide = Chem.MolFromSmiles("CN[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O")
    assert product("meC.[C(4,4)]") == Chem.MolToSmiles(disulfide)
    record = audit_monomer(definition("meC"))
    assert record["activation"]["status"] == "match"
    slots = {
        atom.GetIsotope()
        for atom in definition("meC").GetAtoms()
        if atom.GetAtomicNum() == 0
    }
    assert slots == {1, 2, 4}


def test_completed_terminal_metadata_preserves_restoration_and_connection():
    result = quality_for_monomer(definition("Mpa"))
    assert not any(issue["code"] == "missing_leaving_metadata" for issue in result["issues"])
    assert definition("Mpa").GetProp("m_Rgroups").split(",")[3].strip() == "[H]"
    expected = Chem.MolFromSmiles("OC(=O)[C@H](C)S")
    assert product("Mpa") == Chem.MolToSmiles(expected)
    disulfide = Chem.MolFromSmiles("OC(=O)[C@H](C)SSC[C@H](N)C(=O)O")
    assert product("Mpa.[C(4,4)]") == Chem.MolToSmiles(disulfide)


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


def test_identity_review_is_bound_to_the_reviewed_structure_and_name():
    original = definition("Pca")
    issue = next(i for i in audit_monomer(original)["issues"] if i["code"] == "identity_review")
    assert issue["severity"] == "warning"
    assert "3-carboxylic" in issue["message"]
    renamed = Chem.Mol(original)
    renamed.SetProp("m_name", "Pyridine-3-carboxylate")
    assert not any(i["code"] == "identity_review" for i in audit_monomer(renamed)["issues"])
    assert quality_for_monomer(renamed)["status"] == "unreviewed"


def test_metadata_cannot_claim_a_site_missing_from_the_structure():
    molecule = definition("A")
    molecule.SetProp("m_chem_types", molecule.GetProp("m_chem_types") + ",9:amine_primary")
    issue = next(i for i in audit_monomer(molecule)["issues"] if i["code"] == "orphan_attachment_metadata")
    assert issue["severity"] == "warning"
    assert issue["slots"] == [9]


@pytest.mark.parametrize(
    "symbol,expected_smiles",
    [
        ("Hsl", "N[C@H]1CCOC1=O"),
        ("TATA", "CCC(=O)N1CN(C(=O)CC)CN(C(=O)CC)C1"),
    ],
)
def test_authored_orientation_keeps_the_existing_standalone_product(
    symbol, expected_smiles
):
    result = quality_for_monomer(definition(symbol))
    assert result["activation_status"] == "match"
    assert not any(issue["code"] == "activation_exception" for issue in result["issues"])
    assert product(symbol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected_smiles))


def test_second_hydrogen_role_is_not_a_lysine_chemistry_conflict():
    result = quality_for_monomer(definition("K"))
    assert not any(item["code"] == "chemistry_declaration_difference" for item in result["issues"])
    expected = Chem.MolFromSmiles("N[C@@H](CCCCNC(C)=O)C(=O)O")
    assert product("K.Ac(4,2)") == Chem.MolToSmiles(expected)
    assert product("K.Ac(5,2)") == Chem.MolToSmiles(expected)


@pytest.mark.parametrize("symbol", list("ARNDCEQGHILKMFPSTWYV"))
def test_standard_amino_acids_have_no_library_warnings(symbol):
    result = quality_for_monomer(definition(symbol))
    assert result["status"] == "no_known_exception"
    assert not any(issue["severity"] == "warning" for issue in result["issues"])
