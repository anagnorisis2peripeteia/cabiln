"""Saved work must retain its meaning across library and runtime changes."""

from copy import deepcopy

import pytest
from fastapi.testclient import TestClient

from pyPept.web import projects
from pyPept.web.app import create_app


def document(text="A-G", notation="cabiln"):
    return {"text": text, "notation": notation, "warning": ""}


def project():
    return {
        "format": "cabiln-project",
        "version": 1,
        "document": document(),
        "drafts": {},
        "reference": {"text": "NCC(=O)O", "original": None},
    }


@pytest.fixture
def client():
    with TestClient(create_app(allow_registration=False)) as value:
        yield value


def test_save_open_retains_sources_drafts_and_quality(client):
    value = project()
    converted = client.post("/smiles_to_cabiln", json={"smiles": "NC(C)C(=O)O"}).json()
    value["document"] = {
        **document(converted["cabiln"]),
        "context": converted["context"],
        "quality": {
            key: converted[key]
            for key in (
                "recognition_status",
                "search_complete",
                "inferred_stereo",
                "assignments",
                "warnings",
                "synthetic_components",
            )
        },
    }
    value["drafts"]["cabiln"] = document("K.[unfinished")
    value["reference"]["original"] = {
        "kind": "text",
        "content": "NCC(=O)O",
        "name": "reference.txt",
    }
    response = client.post("/prepare_project", json={"project": value})
    assert response.status_code == 200, response.text
    saved = response.json()["project"]
    assert saved["document"]["quality"] == value["document"]["quality"]
    assert saved["reference"]["text"] == value["reference"]["text"]
    assert saved["reference"]["original"] == value["reference"]["original"]
    assert saved["reference"]["context"] == saved["context"]
    assert saved["drafts"]["cabiln"]["text"] == "K.[unfinished"
    assert saved["context"]["resolutions"]["draft:cabiln"]["resolved"] is None
    assert client.post("/validate_project", json={"project": saved}).json()["valid"]


def test_source_change_is_detected_when_rebinding_library(monkeypatch):
    saved = projects.validate_project(project(), preparing=True)
    current = projects.project_context()
    current["library_binding"]["monomers"] = "new-library"
    monkeypatch.setattr(projects, "project_context", lambda: deepcopy(current))
    # Adding an unrelated definition is allowed if selected chemistry is unchanged.
    assert projects.validate_project(saved)["context"]["library_binding"] == (
        current["library_binding"]
    )
    saved["document"]["text"] = "A-A"
    with pytest.raises(projects.ProjectError, match="definitions differ"):
        projects.validate_project(saved)


@pytest.mark.parametrize("change", ["rules", "convention", "unresolved"])
def test_incompatible_environment_never_silently_rebinds(monkeypatch, change):
    value = project()
    if change == "unresolved":
        value["drafts"]["biln"] = document("unknown-symbol", "biln")
    saved = projects.validate_project(value, preparing=True)
    current = projects.project_context()
    if change == "rules":
        current["library_binding"]["reactions"] = "changed-rules"
    elif change == "convention":
        current["canonical"]["rdkit_version"] = "different-version"
    else:
        current["library_binding"]["monomers"] = "changed-library"
    monkeypatch.setattr(projects, "project_context", lambda: deepcopy(current))
    with pytest.raises(projects.ProjectError) as failure:
        projects.validate_project(saved)
    assert failure.value.code == "project_binding_mismatch"


def test_edited_draft_does_not_reuse_old_signature_after_library_change(monkeypatch):
    value = project()
    value["drafts"]["cabiln"] = document("K-A")
    saved = projects.validate_project(value, preparing=True)
    current = projects.project_context()
    current["library_binding"]["aliases"] = "different-aliases"
    monkeypatch.setattr(projects, "project_context", lambda: deepcopy(current))
    saved["drafts"]["cabiln"]["text"] = "K-G"
    with pytest.raises(projects.ProjectError):
        projects.validate_project(saved, preparing=True)


def test_stale_document_context_cannot_hide_in_new_project_context():
    value = project()
    stale = projects.project_context()
    stale["library_binding"]["monomers"] = "old-library"
    value["document"]["context"] = stale
    with pytest.raises(projects.ProjectError):
        projects.validate_project(value, preparing=True)


def test_library_change_during_prepare_rejected(monkeypatch):
    current = projects.project_context()
    after = deepcopy(current)
    after["library_binding"]["monomers"] = "changed-during-save"
    contexts = iter([current, after])
    monkeypatch.setattr(projects, "project_context", lambda: next(contexts))
    with pytest.raises(projects.ProjectError, match="changed while checking"):
        projects.validate_project(project(), preparing=True)


@pytest.mark.parametrize("mutation", ["version", "format", "depth", "size", "draft"])
def test_invalid_project_rejected_before_replacing_work(client, mutation):
    value = project()
    if mutation == "version":
        value["version"] = True
    elif mutation == "format":
        value["format"] = "unknown"
    elif mutation == "depth":
        node = value
        for _ in range(34):
            node["extra"] = {}
            node = node["extra"]
    elif mutation == "size":
        value["extra"] = "x" * projects.MAX_PROJECT_BYTES
    else:
        value["drafts"]["cabiln"] = document("A", "smiles")
    response = client.post("/prepare_project", json={"project": value})
    assert response.status_code == 409
    assert response.json()["code"] == "project_invalid"


def test_chemistry_metadata_requires_its_original_context():
    value = project()
    value["document"]["canonical"] = {"version": 1}
    with pytest.raises(projects.ProjectError, match="Render this document"):
        projects.validate_project(value, preparing=True)


def test_reference_signature_preserves_original_mol_content():
    from rdkit import Chem

    molecule = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    reference = {
        "text": "NCC(=O)O",
        "original": {"kind": "mol", "content": Chem.MolToMolBlock(molecule)},
    }
    original = projects._reference_resolution(reference)
    reference["text"] = "CCO"
    changed = projects._reference_resolution(reference)
    assert changed["source"] != original["source"]
    assert changed["resolved"] == original["resolved"]
    reference["original"] = None
    assert projects._reference_resolution(reference) != original


@pytest.mark.parametrize("changed", ["document", "draft", "reference"])
def test_tampered_saved_sources_rejected_even_in_same_environment(changed):
    value = project()
    value["drafts"]["cabiln"] = document("G")
    value["reference"]["original"] = {"kind": "text", "content": "NCC(=O)O"}
    saved = projects.validate_project(value, preparing=True)
    if changed == "document":
        saved["document"]["text"] = "W"
    elif changed == "draft":
        saved["drafts"]["cabiln"]["text"] = "W"
    else:
        saved["reference"]["text"] = "N[C@@H](C)C(=O)O"
    with pytest.raises(projects.ProjectError, match="sources differ"):
        projects.validate_project(saved)


def test_new_save_allows_edited_drafts_after_clearing_old_chemistry_evidence():
    value = projects.validate_project(project(), preparing=True)
    value["document"]["text"] = "W"
    value["document"]["canonical"] = {"version": 1}
    with pytest.raises(projects.ProjectError, match="text changed"):
        projects.validate_project(value, preparing=True)
    value["document"]["canonical"] = None
    saved = projects.validate_project(value, preparing=True)
    assert projects.validate_project(saved)["document"]["text"] == "W"


def test_real_library_addition_preserves_project_but_selected_structure_change_does_not(
    tmp_path, monkeypatch
):
    from rdkit import Chem

    from pyPept import monomer_store

    _, shipped = monomer_store._load_sdf()
    glycine, alanine = Chem.Mol(shipped["G"]), Chem.Mol(shipped["A"])
    path = tmp_path / "monomers.sdf"
    path.write_text(Chem.SDWriter.GetText(glycine, kekulize=False))
    monkeypatch.setattr(monomer_store, "_SDF_PATH", path)
    monomer_store._invalidate_sdf()
    try:
        value = project()
        value["document"] = document("G")
        saved = projects.validate_project(value, preparing=True)
        monomer_store.register_molecule(alanine)
        assert projects.validate_project(saved)["document"]["text"] == "G"
        substituted = Chem.Mol(alanine)
        substituted.SetProp("m_abbr", "G")
        substituted.SetProp("symbol", "G")
        path.write_text(
            Chem.SDWriter.GetText(substituted, kekulize=False)
            + Chem.SDWriter.GetText(alanine, kekulize=False)
        )
        monomer_store._invalidate_sdf()
        with pytest.raises(projects.ProjectError, match="definitions differ"):
            projects.validate_project(saved)
    finally:
        monomer_store._invalidate_sdf()


def test_reference_keeps_independent_binding_when_main_context_changes(monkeypatch):
    saved = projects.validate_project(project(), preparing=True)
    current = projects.project_context()
    current["library_binding"]["monomers"] = "new-library"
    monkeypatch.setattr(projects, "project_context", lambda: deepcopy(current))
    saved["context"] = deepcopy(current)
    saved["document"]["context"] = deepcopy(current)
    # The old reference proof still demonstrates unchanged reference chemistry.
    prepared = projects.validate_project(saved, preparing=True)
    assert prepared["reference"]["context"]["library_binding"] == (
        current["library_binding"]
    )
    resolve = projects._reference_resolution

    def changed(reference):
        signature = resolve(reference)
        signature["resolved"] = "changed-reference-chemistry"
        return signature

    monkeypatch.setattr(projects, "_reference_resolution", changed)
    with pytest.raises(projects.ProjectError, match="definitions differ"):
        projects.validate_project(saved, preparing=True)


@pytest.mark.parametrize(
    "route,payload",
    [
        ("/render_reference", {"input": "G-G"}),
        ("/render_smiles", {"smiles": "NCC(=O)O"}),
    ],
)
def test_reference_render_returns_actual_library_context(client, route, payload):
    rendered = client.post(route, json=payload)
    assert rendered.status_code == 200, rendered.text
    assert rendered.json()["context"] == projects.project_context()
