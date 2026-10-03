"""Ingestion paths preserve the same usable, persistent monomer definitions."""

import csv

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import monomer_store
from pyPept.interfaces.cli_monomer import register_monomer
from pyPept.interfaces.monomer_pipeline import (
    build_library_from_csv,
    pre_activate,
    write_sdf,
)
from pyPept.web.app import create_app


def author_csv(path, row):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def one_record(path):
    records = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert len(records) == 1
    assert records[0] is not None
    return records[0]


@pytest.mark.parametrize("source", ["NCC(=O)O", "N[C@@H](Cc1ccc(O)cc1)C(=O)O"])
def test_unchanged_csv_rebuild_preserves_derived_chemistry(tmp_path, source):
    csv_path, sdf_path = tmp_path / "monomers.csv", tmp_path / "monomers.sdf"
    author_csv(csv_path, {"token": "Imported", "input": source, "type": "aa"})
    assert build_library_from_csv(csv_path, sdf_path) == (1, [])
    original = one_record(sdf_path)
    declared = original.GetProp("m_chem_types")
    assert monomer_store.parse_chem_types(declared) == pre_activate(source).chem_types

    assert build_library_from_csv(csv_path, sdf_path) == (1, [])
    rebuilt = one_record(sdf_path)
    assert rebuilt.GetProp("m_chem_types") == declared
    assert rebuilt.GetProp("m_Rgroups") == original.GetProp("m_Rgroups")
    assert rebuilt.GetProp("m_activation_policy") == original.GetProp("m_activation_policy") == 'canonical-sites-v1'
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(original)


def test_csv_rebuild_preserves_nonstandard_slots_and_author_metadata(tmp_path):
    csv_path, sdf_path = tmp_path / "monomers.csv", tmp_path / "monomers.sdf"
    template = "[1*]N([3*])C(CO[7*])C([2*])=O"
    declaration = "1:backbone_n,2:backbone_c,3:backbone_n_mod,7:hydroxyl"
    author_csv(
        csv_path,
        {
            "token": "ExtraSlot",
            "input": template,
            "chuckles": template,
            "type": "aa",
            "r1_leaving": "[H]",
            "r2_leaving": "[OH]",
            "r3_leaving": "[H]",
            "r7_leaving": "[H]",
            "chem_types": declaration,
            "author_note": "Keep this experimental annotation",
        },
    )
    for _ in range(2):
        assert build_library_from_csv(csv_path, sdf_path) == (1, [])
        stored = one_record(sdf_path)
        assert stored.GetProp("m_Rgroups").split(",")[6] == "[H]"
        assert stored.GetProp("m_chem_types") == declaration
        with csv_path.open(newline="") as stream:
            row = next(csv.DictReader(stream))
        assert row["r7_leaving"] == "[H]"
        assert row["author_note"] == "Keep this experimental annotation"


def test_bulk_writer_preserves_sparse_declarations_and_standalone_records(tmp_path):
    path = tmp_path / "monomers.sdf"
    rows = [
        {
            "token": "Sparse",
            "chuckles": "[1*]NCC([2*])=O",
            "r1_leaving": "[H]",
            "r2_leaving": "[OH]",
            "_chem_types": {2: "backbone_c"},
        },
        {"token": "Standalone", "chuckles": "O"},
    ]
    assert write_sdf(rows, path) == (2, [])
    sparse, standalone = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert sparse.GetProp("m_chem_types") == "2:backbone_c"
    assert Chem.MolToSmiles(standalone) == "O"
    assert standalone.GetProp("m_Rgroups") == "None,None,None"


def test_bulk_writer_preserves_legacy_incomplete_slot_metadata(tmp_path):
    path = tmp_path / "monomers.sdf"
    rows = [
        {"token": "LegacyCap", "chuckles": "[2*]C(=O)OC(C)(C)C", "r1_leaving": "[Cl]"},
        {"token": "LegacyThiol", "chuckles": "[2*]C(=O)CS[7*]", "r2_leaving": "[OH]"},
    ]
    assert write_sdf(rows, path) == (2, [])
    cap, thiol = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert cap.GetProp("m_Rgroups") == "[Cl]"
    assert thiol.GetProp("m_Rgroups") == "None,[OH]"


def test_bulk_build_reports_record_validation_errors(tmp_path):
    csv_path, sdf_path = tmp_path / "monomers.csv", tmp_path / "monomers.sdf"
    author_csv(
        csv_path,
        {"token": "InvalidLeaving", "chuckles": "[1*]NCC([2*])=O", "r1_leaving": "CO"},
    )
    count, errors = build_library_from_csv(csv_path, sdf_path)
    assert count == 0
    assert len(errors) == 1
    assert "InvalidLeaving" in errors[0] and "single atom" in errors[0]


def test_record_construction_retains_unrelated_properties_without_mutation():
    activated = pre_activate("NCC(=O)O")
    original = Chem.MolFromSmiles(activated.chuckles)
    original.SetProp("author_note", "Independent library annotation")
    stored = monomer_store.monomer_record(
        original, "Annotated", activated.leaving, activated.chem_types
    )
    assert stored.GetProp("author_note") == original.GetProp("author_note")
    assert not original.HasProp("symbol")
    assert stored.GetProp("symbol") == "Annotated"


def test_new_registration_rejects_contradictory_chemistry_before_writing(tmp_path, monkeypatch):
    path = tmp_path / 'monomers.sdf'
    path.touch()
    monkeypatch.setenv('CABILN_MONOMER_LIBRARY', str(path))
    activated = pre_activate('N[C@@H](CS)C(=O)O')
    with TestClient(create_app(allow_registration=True)) as client:
        response = client.post('/register_monomer', json={
            'abbr': 'FalseThiol', 'name': 'Contradictory declaration',
            'chuckles': activated.chuckles, 'leaving': activated.leaving,
            'chem_types': {**activated.chem_types, 4: 'hydroxyl'},
        })
    assert response.status_code == 400, response.text
    assert 'R4 is thiol' in response.text
    assert path.read_bytes() == b''


@pytest.mark.parametrize("entry", ["csv", "web"])
def test_extended_smiles_keeps_numbered_slots_through_registration(
    tmp_path, monkeypatch, entry
):
    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence

    path = tmp_path / "monomers.sdf"
    path.touch()
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    template = "[1*]N[C@@H]([13CH3])C([2*])=O |a:2|"
    try:
        if entry == "csv":
            csv_path = tmp_path / "monomers.csv"
            author_csv(
                csv_path,
                {
                    "token": "AbsoluteAla",
                    "input": template,
                    "type": "aa",
                    "r1_leaving": "[H]",
                    "r2_leaving": "[OH]",
                },
            )
            for _ in range(2):
                assert build_library_from_csv(csv_path, path) == (1, [])
        else:
            with TestClient(create_app(allow_registration=True)) as client:
                response = client.post(
                    "/register_monomer",
                    json={
                        "abbr": "AbsoluteAla",
                        "name": "Absolute alanine",
                        "chuckles": template,
                        "leaving": {1: "[H]", 2: "[OH]"},
                        "chem_types": {1: "backbone_n", 2: "backbone_c"},
                    },
                )
                assert response.status_code == 200, response.text
        stored = one_record(path)
        assert {
            atom.GetIsotope() for atom in stored.GetAtoms() if atom.GetAtomicNum() == 0
        } == {1, 2}
        assert stored.GetProp("m_chem_types") == "1:backbone_n,2:backbone_c"
        product = Molecule(Sequence("AbsoluteAla")).mol
        assert Chem.MolToSmiles(product) == Chem.MolToSmiles(
            Chem.MolFromSmiles("N[C@@H]([13CH3])C(=O)O")
        )
    finally:
        monomer_store._invalidate_sdf()


def test_cli_and_http_construct_equivalent_definitions(tmp_path, monkeypatch):
    path = tmp_path / "monomers.sdf"
    path.touch()
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    source = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"
    try:
        register_monomer(source, "CliMonomer", name="Same definition")
        with TestClient(create_app(allow_registration=True)) as client:
            preview = client.post("/preview_monomer", json={"smiles": source})
            assert preview.status_code == 200, preview.text
            activated = preview.json()
            response = client.post(
                "/register_monomer",
                json={
                    **{
                        key: activated[key]
                        for key in ("chuckles", "chem_types", "leaving")
                    },
                    "abbr": "WebMonomer",
                    "name": "Same definition",
                    "type": "aa",
                    "subtype": "modified",
                },
            )
            assert response.status_code == 200, response.text
        cli, web = list(Chem.SDMolSupplier(str(path), removeHs=False))
        assert Chem.MolToSmiles(cli) == Chem.MolToSmiles(web)
        for field in ("m_name", "m_type", "m_subtype", "m_chem_types"):
            assert cli.GetProp(field) == web.GetProp(field)
        # Each adapter retains its historical padding convention.
        cli_groups = cli.GetProp("m_Rgroups").split(",")
        web_groups = web.GetProp("m_Rgroups").split(",")
        assert web_groups[: len(cli_groups)] == cli_groups
        assert all(group == "None" for group in web_groups[len(cli_groups) :])
    finally:
        monomer_store._invalidate_sdf()
