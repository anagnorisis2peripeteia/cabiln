"""Public chemistry boundaries must retain stereo and selectable definitions."""

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from rdkit import Chem

from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.web.app import create_app
from pyPept.web.schemas import _RegisterReq


@pytest.fixture
def client():
    with TestClient(create_app(allow_registration=False)) as client:
        yield client


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_relative_stereo_is_rejected_instead_of_claiming_a_defined_conversion(
    client, group
):
    source = f"N[C@@H](C)C(=O)N[C@@H](C)C(=O)O |{group}:1,6|"
    response = client.post("/smiles_to_cabiln", json={"smiles": source})
    assert response.status_code == 400
    assert "AND/OR" in response.json()["error"]


@pytest.mark.parametrize("group", ["&1", "o1"])
@pytest.mark.parametrize("route", ["/render_reference", "/render_mol"])
def test_normalized_reference_cannot_discard_unsupported_stereo(client, group, route):
    source = f"N[C@@H](C)C(=O)N[C@@H](C)C(=O)O |{group}:1,6|"
    molecule = Chem.MolFromSmiles(source)
    assert len(molecule.GetStereoGroups()) == 1
    payload = (
        {"input": source}
        if route == "/render_reference"
        else {"mol_block": Chem.MolToMolBlock(molecule, forceV3000=True)}
    )
    response = client.post(route, json=payload)
    assert response.status_code == 400, response.text
    assert "AND/OR" in response.json()["error"]
    assert "smiles" not in response.json()


@pytest.mark.parametrize("suffix", ["", " |a:1,6|"])
@pytest.mark.parametrize("route", ["/render_reference", "/render_mol"])
def test_supported_reference_normalization_preserves_exact_stereo(
    client, suffix, route
):
    source = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)O" + suffix
    molecule = Chem.MolFromSmiles(source)
    payload = (
        {"input": source}
        if route == "/render_reference"
        else {"mol_block": Chem.MolToMolBlock(molecule, forceV3000=True)}
    )
    response = client.post(route, json=payload)
    assert response.status_code == 200, response.text
    normalized = response.json()["smiles"]
    assert normalized == Chem.MolToSmiles(molecule)
    converted = client.post("/smiles_to_cabiln", json={"smiles": normalized})
    assert converted.status_code == 200, converted.text
    assert converted.json()["inferred_stereo"] is False
    verified = client.post(
        "/verify", json={"smiles": source, "cabiln": converted.json()["cabiln"]}
    )
    assert verified.status_code == 200, verified.text
    assert verified.json()["match"] is True


def test_every_displayed_symbol_passes_its_http_selection_boundary(client, monkeypatch):
    from pyPept.editor import EditResult
    from pyPept.web import builder

    class SelectionBoundary:
        """Let HTTP validation run without asserting arbitrary reaction pairs."""

        def __init__(self, source):
            pass

        def select(self, index):
            return index

        def insert_backbone(self, selected, symbol):
            return symbol

        def attach(self, host, r_host, symbol, r_new, *, with_details=False):
            return EditResult(symbol, (), (), ()) if with_details else symbol

    response = client.get("/monomers")
    assert response.status_code == 200, response.text
    entries = response.json()
    assert "L_hArg(Et,Et)" in {entry["abbr"] for entry in entries}
    monkeypatch.setattr(builder, "PeptideDocument", SelectionBoundary)
    for entry in entries:
        if entry["backbone_insertable"]:
            selected = client.post(
                "/insert_backbone",
                json={"cabiln": "A-G", "after_idx": 0, "new_abbr": entry["abbr"]},
            )
            assert selected.status_code == 200, (entry["abbr"], selected.text)
            assert selected.json()["result"] == entry["abbr"]
        for symbol in {
            entry[key] for key in ("abbr", "nterm_abbr", "cterm_abbr") if key in entry
        }:
            selected = client.post(
                "/insert_bond",
                json={
                    "cabiln": "K-A",
                    "host_residue_idx": 0,
                    "new_abbr": symbol,
                    "r_host": 4,
                    "r_new": 2,
                },
            )
            assert selected.status_code == 200, (symbol, selected.text)
            assert selected.json()["result"] == symbol


@pytest.mark.parametrize(
    "route,payload,expected",
    [
        (
            "/insert_backbone",
            {"cabiln": "", "after_idx": 0},
            "L_hArg(Et,Et)",
        ),
        (
            "/insert_backbone",
            {"cabiln": "A-G", "after_idx": 0},
            "A-L_hArg(Et,Et)-G",
        ),
        (
            "/insert_bond",
            {"cabiln": "K-A", "host_residue_idx": 0, "r_host": 4, "r_new": 2},
            "K.{L_hArg(Et,Et)(4,2)}-A",
        ),
    ],
)
def test_existing_punctuated_definition_is_insertable(client, route, payload, expected):
    response = client.post(route, json={**payload, "new_abbr": "L_hArg(Et,Et)"})
    assert response.status_code == 200, response.text
    assert response.json()["result"] == expected
    rendered = client.post("/render", json={"cabiln": expected})
    assert rendered.status_code == 200, rendered.text


def test_registration_names_remain_stricter_than_existing_selection():
    with pytest.raises(ValidationError, match="abbr"):
        _RegisterReq(
            chuckles="[1*]NCC([2*])=O",
            chem_types={},
            leaving={},
            abbr="New(Et,Et)",
            name="A new record",
        )


def test_high_numbered_site_can_be_inspected_and_connected(client):
    source = "<[1*]N[C@@H](CS[65*])C([2*])=O>"
    inspected = client.get(
        "/monomer_rgroups", params={"abbr": "", "residue_idx": 0, "cabiln": source}
    )
    assert inspected.status_code == 200, inspected.text
    assert 65 in {site["slot"] for site in inspected.json()["rgroups"]}
    payload = {
        "cabiln": source,
        "host_residue_idx": 0,
        "new_abbr": "C",
        "r_host": 65,
        "r_new": 4,
    }
    response = client.post("/insert_bond", json=payload)
    assert response.status_code == 200, response.text
    product = Molecule(Sequence(response.json()["result"]), depiction=None).mol
    expected = Chem.MolFromSmiles("N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O")
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(expected)
    invalid = client.post("/insert_bond", json={**payload, "r_host": 66})
    assert invalid.status_code == 400
    assert "R66" in invalid.json()["error"]


def test_selection_still_requires_one_resolved_monomer(client):
    for invalid in ("MissingLaunchMonomer", "A-G"):
        response = client.post(
            "/insert_backbone",
            json={"cabiln": "A-G", "after_idx": 0, "new_abbr": invalid},
        )
        assert response.status_code == 400, response.text


def test_render_context_and_residue_kind_follow_the_resolved_document(client):
    from pyPept.canonical import canonical_convention
    from pyPept.monomer_store import library_binding

    payload = {"cabiln": "A-<[1*]NCC([2*])=O>"}
    expected = {
        "project_version": 1,
        "library_binding": library_binding(),
        "canonical": canonical_convention(),
    }
    for _ in range(2):
        response = client.post("/render", json=payload)
        assert response.status_code == 200, response.text
        assert response.json()["context"] == expected
        assert [residue["kind"] for residue in response.json()["residues"]] == [
            "library",
            "synthetic",
        ]


def test_render_rejects_a_library_change_before_returning_success(client, monkeypatch):
    from pyPept import monomer_store
    from pyPept.web import rendering

    binding = monomer_store.library_binding()
    changed = {**binding, "monomers": "changed"}
    calls = iter((binding, changed))
    monkeypatch.setattr(monomer_store, "library_binding", lambda: next(calls))
    monkeypatch.setattr(rendering, "_rc_get", lambda _: None)
    response = client.post("/render", json={"cabiln": "A-G"})
    assert response.status_code == 400, response.text
    assert "library changed" in response.json()["error"]


@pytest.mark.parametrize(
    "route,payload",
    [
        ("/render", {"cabiln": "A-G"}),
        ("/render_smiles", {"smiles": "NCC(=O)O"}),
        ("/render_reference", {"input": "NCC(=O)O"}),
        (
            "/render_mol",
            {"mol_block": Chem.MolToMolBlock(Chem.MolFromSmiles("NCC(=O)O"))},
        ),
    ],
)
def test_render_internal_failures_are_distinct_and_sanitized(
    client, monkeypatch, route, payload
):
    from pyPept.web import rendering

    def fail(*args, **kwargs):
        raise RuntimeError("private submitted structure and implementation details")

    monkeypatch.setattr(rendering, "_rc_get", lambda _: None)
    monkeypatch.setattr(rendering, "_draw_mol", fail)
    response = client.post(route, json=payload)
    assert response.status_code == 500, response.text
    assert response.json()["error"] == "An internal error occurred."
    assert "private" not in response.text
    assert response.json()["request_id"]
