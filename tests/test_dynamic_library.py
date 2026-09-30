"""An unseen monomer must flow from detection through palette to assembly."""

import pytest
from hypothesis import given, strategies as st
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import monomer_store
from pyPept.interfaces.cli_monomer import register_monomer
from pyPept.web.app import create_app

from _chemistry_fuzz import alternate_smiles, isolated_library, reordered, write_library
from _fuzzing import fuzz_settings, record


@pytest.fixture
def library_app(tmp_path, monkeypatch):
    _, existing = monomer_store._load_sdf()
    path = tmp_path / "monomers.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for symbol in ("ac", "am", "G", "C"):
            writer.write(existing[symbol])
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    with TestClient(create_app(allow_registration=True)) as client:
        yield path, client
    monomer_store._invalidate_sdf()


def post(client, route, payload):
    response = client.post(route, json=payload)
    assert response.status_code == 200, response.text
    return response.json()


def smiles_of_render(client, notation):
    data = post(client, "/render", {"cabiln": notation})
    return Chem.MolToSmiles(Chem.MolFromMolBlock(data["mol_block"]))


def test_new_monomer_is_detected_listed_selected_and_bonded(library_app):
    path, client = library_app
    # Warm all relevant caches before ingesting a previously unseen symbol.
    assert len(client.get("/monomers").json()) == 4
    smiles_of_render(client, "ac-G-am")
    source = "N[C@@H](CCCS)C(=O)O"
    preview = post(client, "/preview_monomer", {"smiles": source})
    assert preview["chem_types"]["4"] == "thiol"
    post(
        client,
        "/register_monomer",
        {
            **{key: preview[key] for key in ("chuckles", "chem_types", "leaving")},
            "abbr": "NewThiol",
            "name": "Unseen thiol monomer",
        },
    )
    palette = client.get("/monomers").json()
    tile = next(item for item in palette if item["abbr"] == "NewThiol")
    assert "4:thiol" in tile["chem_types"]
    assert client.get("/monomer_svg", params={"abbr": "NewThiol"}).status_code == 200
    assert smiles_of_render(client, "NewThiol") == Chem.MolToSmiles(
        Chem.MolFromSmiles(source)
    )
    inserted = post(
        client,
        "/insert_backbone",
        {
            "cabiln": "ac-G-am",
            "after_idx": 1,
            "new_abbr": "NewThiol",
        },
    )["result"]
    assert smiles_of_render(client, inserted) == Chem.MolToSmiles(
        Chem.MolFromSmiles("CC(=O)NCC(=O)N[C@@H](CCCS)C(=O)N")
    )
    slots = client.get(
        "/monomer_rgroups",
        params={
            "abbr": "NewThiol",
            "cabiln": inserted,
            "residue_idx": 2,
        },
    ).json()["rgroups"]
    assert {site["slot"] for site in slots if site["used"]} == {1, 2}
    assert next(site for site in slots if site["slot"] == 4)["chem_type"] == "thiol"
    connected = post(
        client,
        "/insert_bond",
        {
            "cabiln": inserted,
            "host_residue_idx": 2,
            "new_abbr": "C",
            "r_host": 4,
            "r_new": 4,
        },
    )["result"]
    assert smiles_of_render(client, connected) == Chem.MolToSmiles(
        Chem.MolFromSmiles("CC(=O)NCC(=O)N[C@@H](CCCSSC[C@H](N)C(=O)O)C(=O)N")
    )
    assert len(list(Chem.SDMolSupplier(str(path)))) == 5


def test_cli_ingestion_updates_running_palette_and_assembly(library_app):
    _, client = library_app
    initial = client.get("/monomers").json()
    source = "N[C@@H](CF)C(=O)O"
    register_monomer(source, "NewFluoro")
    updated = client.get("/monomers").json()
    assert len(updated) == len(initial) + 1
    assert "NewFluoro" in {item["abbr"] for item in updated}
    assert smiles_of_render(client, "NewFluoro") == Chem.MolToSmiles(
        Chem.MolFromSmiles(source)
    )


def test_detection_only_labels_from_preview_can_be_registered(library_app):
    _, client = library_app
    preview = post(
        client, "/preview_monomer", {"smiles": "N[C@@H](Cc1ccc(O)cc1)C(=O)O"}
    )
    assert preview["chem_types"]["4"] == "hydroxyl_phenolic"
    post(
        client,
        "/register_monomer",
        {
            **{key: preview[key] for key in ("chuckles", "chem_types", "leaving")},
            "abbr": "NewPhenol",
            "name": "Detected phenol",
        },
    )
    tile = next(
        item for item in client.get("/monomers").json() if item["abbr"] == "NewPhenol"
    )
    slots = client.get("/monomer_rgroups", params={"abbr": "NewPhenol"}).json()[
        "rgroups"
    ]
    effective = ",".join(f"{site['slot']}:{site['chem_type']}" for site in slots)
    assert tile["chem_types"] == effective


def test_sdf_record_without_declared_types_is_discovered_and_detected(library_app):
    path, client = library_app
    client.get("/monomers")
    register_monomer("N[C@@H](CCCS)C(=O)O", "ImportedThiol")
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    molecules[-1].ClearProp("m_chem_types")
    with Chem.SDWriter(str(path)) as writer:
        for molecule in molecules:
            writer.write(molecule)
    tile = next(
        item
        for item in client.get("/monomers").json()
        if item["abbr"] == "ImportedThiol"
    )
    assert "4:thiol" in tile["chem_types"]
    slots = client.get("/monomer_rgroups", params={"abbr": "ImportedThiol"}).json()[
        "rgroups"
    ]
    assert next(site for site in slots if site["slot"] == 4)["chem_type"] == "thiol"


@pytest.mark.fuzz
@fuzz_settings(examples=10)
@given(length=st.integers(2, 4), stereo=st.sampled_from(("@", "@@")),
       isotope=st.sampled_from(("", "13")), order=st.integers(0, 65535))
def test_fuzz_registration_refreshes_consumers_and_preserves_selected_bindings(
    length, stereo, isotope, order
):
    from pyPept.interfaces.monomer_pipeline import pre_activate
    from pyPept.smiles import convert_smiles
    from pyPept.web import projects

    source = f"N[{isotope}C{stereo}H]({'C' * length}S)C(=O)O"
    record("library.register-reorder-edit", source=source, order=order)

    def save(text):
        return projects.validate_project({
            "format": "cabiln-project", "version": 1,
            "document": {"text": text, "notation": "cabiln", "warning": ""},
            "drafts": {}, "reference": {"text": "NCC(=O)O", "original": None},
        }, preparing=True)

    with isolated_library(("ac", "am", "G", "C")) as path:
        with TestClient(create_app(allow_registration=True)) as client:
            assert len(client.get("/monomers").json()) == 4
            baseline = smiles_of_render(client, "ac-G-am")
            assert convert_smiles("NCC(=O)O").recognition_status == "complete"
            saved = save("ac-G-am")
            preview = post(client, "/preview_monomer", {"smiles": alternate_smiles(source, order)})
            assert preview["chem_types"]["4"] == "thiol"
            post(client, "/register_monomer", {
                **{key: preview[key] for key in ("chuckles", "chem_types", "leaving")},
                "abbr": "Novel", "name": "Generated temporary thiol",
            })
            tile = next(item for item in client.get("/monomers").json() if item["abbr"] == "Novel")
            assert "4:thiol" in tile["chem_types"]
            slots = client.get("/monomer_rgroups", params={"abbr": "Novel"}).json()["rgroups"]
            assert {site["slot"] for site in slots} == {1, 2, 3, 4}
            assert next(site for site in slots if site["slot"] == 4)["chem_type"] == "thiol"
            assert smiles_of_render(client, "Novel") == Chem.MolToSmiles(Chem.MolFromSmiles(source))
            # Product expectation is independent of activation and assembly.
            product = f"CC(=O)NCC(=O)N[{isotope}C{stereo}H]({'C' * length}S)C(=O)N"
            assert smiles_of_render(client, "ac-G-Novel-am") == Chem.MolToSmiles(Chem.MolFromSmiles(product))
            recognized = convert_smiles(alternate_smiles(source, order))
            assert recognized.recognition_status == "complete"
            assert len(recognized.assignments) == 1 and recognized.assignments[0].symbol == "Novel"
            assert projects.validate_project(saved)["document"]["text"] == "ac-G-am"
            selected = save("Novel")
            records = list(Chem.SDMolSupplier(str(path), removeHs=False))
            # Changing physical order and atom order must not change chemistry.
            write_library(path, (reordered(molecule, order) for molecule in reversed(records)))
            assert smiles_of_render(client, "ac-G-am") == baseline
            assert projects.validate_project(selected)["document"]["text"] == "Novel"
            assert convert_smiles(source).recognition_status == "complete"

            edited_source = f"N[{isotope}C{stereo}H]({'C' * length}O)C(=O)O"
            activated = pre_activate(edited_source)
            replacement = monomer_store.monomer_record(
                activated.chuckles, "Novel", activated.leaving, activated.chem_types)
            write_library(path, (replacement if item.GetProp("m_abbr") == "Novel" else item for item in records))
            with pytest.raises(projects.ProjectError, match="definitions differ"):
                projects.validate_project(selected)
            assert projects.validate_project(saved)["document"]["text"] == "ac-G-am"
            assert smiles_of_render(client, "Novel") == Chem.MolToSmiles(Chem.MolFromSmiles(edited_source))
            new_tile = next(item for item in client.get("/monomers").json() if item["abbr"] == "Novel")
            assert "4:hydroxyl" in new_tile["chem_types"] and "thiol" not in new_tile["chem_types"]
            updated = convert_smiles(edited_source)
            assert updated.recognition_status == "complete"
            assert updated.assignments[0].symbol == "Novel"
            old = convert_smiles(source)
            assert not any(item.recognized and item.symbol == "Novel" for item in old.assignments)
