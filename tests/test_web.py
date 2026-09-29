"""Exercise builder edits, rendering, and library writes through HTTP."""

from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import monomer_store
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.web import rendering
from pyPept.web.app import create_app


@pytest.fixture
def client():
    with TestClient(create_app(allow_registration=False)) as client:
        yield client


def canonical(notation):
    return Chem.MolToSmiles(Molecule(Sequence(notation)).get_molecule(fmt="ROMol"))


@pytest.mark.parametrize(
    "source",
    [
        "ac-K.[[K(4,2).A(1,2)garbage].ac(4,2)]-am",
        "ac-K.[!x(4,1)]-D.!x(4,4)-am",
    ],
)
@pytest.mark.parametrize("target", [None, "branch", "bracket"])
def test_notation_routes_reject_discarded_text_and_conflicting_slots(
    client, source, target
):
    path = "/convert_notation" if target else "/render"
    body = {"cabiln": source}
    if target:
        body["target"] = target
    response = client.post(path, json=body)
    assert response.status_code == 400, response.text
    assert response.json()["error"]
    assert "svg" not in response.json()
    assert "result" not in response.json()
    assert client.get("/health").status_code == 200


@pytest.mark.parametrize("notation", ["percent", "bracket"])
def test_conversion_assignments_follow_rendered_occurrences(client, notation):
    source = canonical("ac-K.!1(4,2)-A-am%ac-G-G.!1(2,4)")
    response = client.post(
        "/smiles_to_cabiln", json={"smiles": source, "notation": notation}
    )
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["recognition_status"] == "complete"
    rendered = client.post("/render", json={"cabiln": data["cabiln"]}).json()
    expected_names = {r["idx"]: r["abbr"] for r in rendered["residues"]}
    assert {
        a["residue_index"]: a["symbol"] for a in data["assignments"]
    } == expected_names
    molecule = Chem.MolFromMolBlock(rendered["mol_block"])
    original = Chem.MolFromSmiles(source)
    assert Chem.MolToSmiles(molecule) == Chem.MolToSmiles(original)
    matches = molecule.GetSubstructMatches(original, useChirality=True, uniquify=False)
    assert any(
        all(
            {match[index] for index in assignment["source_atoms"]}
            == set(rendered["residue_map"][str(assignment["residue_index"])])
            for assignment in data["assignments"]
        )
        for match in matches
    )


def test_unresolved_component_does_not_shift_known_assignment_index(client):
    source = "O=C1NCCCC(C(F)(F)F)C1.NCC(=O)O"
    response = client.post("/smiles_to_cabiln", json={"smiles": source})
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["recognition_status"] == "unresolved"
    assert len(data["assignments"]) == 1
    assignment = data["assignments"][0]
    assert assignment["symbol"] == "G"
    assert assignment["residue_index"] == 1
    assert (
        Sequence(data["cabiln"]).s_monomers[assignment["residue_index"]]["m_abbr"]
        == "G"
    )


@pytest.mark.parametrize("notation", ["A-G", "!1-A-G-A-!1", "K.{G(4,2).ac(1,2)}-A"])
def test_render_exports_same_molecule_and_residue_map(client, notation):
    response = client.post("/render", json={"cabiln": notation})
    assert response.status_code == 200, response.text
    data = response.json()
    mol = Chem.MolFromMolBlock(data["mol_block"])
    assert Chem.MolToSmiles(mol) == canonical(notation)
    assert "<svg" in data["svg"]
    assert data["cabiln_echo"] == notation
    assert set(data["residue_map"]) == {str(r["idx"]) for r in data["residues"]}


@pytest.mark.parametrize(
    "source,roots,groups",
    [
        ("A%K.{G(4,2)}-A", [[0], [1, 2]], [(1, [3], "{", None)]),
        (
            "ac-K.[G(4,2)]-K.{A(4,2)}-am",
            [[0, 1, 2, 3]],
            [(1, [4], "[", None), (2, [5], "{", None)],
        ),
        (
            "ac-K.{G(4,2)[.A(1,2)]}-D-am",
            [[0, 1, 2, 3]],
            [(1, [4], "{", None), (4, [5], "[", 0)],
        ),
    ],
)
def test_render_layout_retains_each_source_group(client, source, roots, groups):
    response = client.post("/render", json={"cabiln": source})
    assert response.status_code == 200, response.text
    data = response.json()
    layout = data["layout"]
    assert [s["roots"] for s in layout["segments"]] == roots
    assert [
        (g["host"], g["roots"], g["opening"], g["parent"]) for g in layout["groups"]
    ] == groups
    displayed = [r for segment in roots for r in segment]
    displayed.extend(r for group in layout["groups"] for r in group["roots"])
    assert sorted(displayed) == list(range(len(data["residues"])))
    assert len(displayed) == len(set(displayed))


@pytest.mark.parametrize(
    "source,expected",
    [
        (
            "!r-C-A-C-!r",
            [(0, "!r", [0, 2], None, True), (2, "!r", [0, 2], None, False)],
        ),
        (
            "C.!r(1,2).!s(4,4)-A-C.{!s(4,4).!r(2,1)}",
            [
                (0, "!r", [0, 2], None, False),
                (0, "!s", [0, 2], None, False),
                (2, "!s", [0, 2], 0, False),
                (2, "!r", [0, 2], 0, False),
            ],
        ),
        (
            "ac-K.{G(4,2)[.A(1,2).!1(1,4)]}-D.!1(4,1)-am",
            [(5, "!1", [2, 5], 1, False), (2, "!1", [2, 5], None, False)],
        ),
        (
            "ac-K.[G(4,2)[.!1(1,4)]]-D.!1(4,1)-am",
            [(4, "!1", [2, 4], 1, False), (2, "!1", [2, 4], None, False)],
        ),
    ],
)
def test_render_markers_refer_to_their_source_owner_and_connection(
    client, source, expected
):
    response = client.post("/render", json={"cabiln": source})
    assert response.status_code == 200, response.text
    assert [
        (m["residue"], m["tag"], m["members"], m["group"], m["before"])
        for m in response.json()["layout"]["markers"]
    ] == expected


def test_legacy_render_returns_the_source_that_its_residue_ids_describe(client):
    raw = "D.(4,1)-G-am%G-A-am%K-A"
    response = client.post("/render", json={"cabiln": raw})
    assert response.status_code == 200, response.text
    data = response.json()
    normalized = data["normalized_cabiln"]
    assert data["cabiln_echo"] == normalized
    assert normalized != raw
    expected = "D.[G(4,1).A(2,1).am(2,1)]-G-am%K-A"
    assert canonical(normalized) == canonical(expected)
    assert [r["abbr"] for r in data["residues"]] == [
        m["m_abbr"] for m in Sequence(normalized).s_monomers
    ]
    assert any("legacy positional" in message for message in data["warnings"])
    conversion = client.post(
        "/convert_notation", json={"cabiln": raw, "target": "bracket"}
    )
    assert conversion.status_code == 200, conversion.text
    assert canonical(conversion.json()["result"]) == canonical(expected)


@pytest.mark.parametrize(
    "route,payload",
    [
        ("/render", {"cabiln": "bogus"}),
        ("/verify", {"cabiln": "bogus", "smiles": "CC"}),
        ("/render_reference", {"input": "bogus"}),
        ("/to_cabiln", {"input": "bogus"}),
        ("/smiles_to_cabiln", {"smiles": "bogus"}),
    ],
)
def test_invalid_chemistry_returns_error_without_stopping_server(
    client, route, payload
):
    response = client.post(route, json=payload)
    assert response.status_code == 400
    assert response.json()["error"]
    assert client.get("/health").json() == {"status": "ok"}


@pytest.mark.parametrize(
    "payload",
    [
        {"cabiln": ""},
        {"cabiln": " "},
        {"cabiln": "A", "width": 1000000000},
        {"cabiln": "A", "height": -1},
        {"cabiln": "A", "seed": -1},
        {"cabiln": "A" * 20001},
    ],
)
def test_request_limits_produce_ui_readable_errors(client, payload):
    response = client.post("/render", json=payload)
    assert response.status_code == 422
    assert response.json()["error"]


def test_query_dimension_limits(client):
    assert client.get("/monomer_svg?abbr=A&width=99999999").status_code == 422


def test_crosslink_renumbering_preserves_ten_and_one(client):
    notation = "C.!2(4,4)-C.!1(4,4)-C.!10(4,4)-C.!10-C.!1-C.!2"
    response = client.post(
        "/convert_notation", json={"cabiln": notation, "target": "bracket"}
    )
    assert response.status_code == 200, response.text
    result = response.json()["result"]
    assert "!20" not in result
    assert canonical(result) == canonical(notation)


def test_notation_conversion_refuses_an_assemblable_but_different_product(
    client, monkeypatch
):
    from dataclasses import replace

    from pyPept.peptide import serialize

    def incorrect(peptide, notation, **policy):
        return replace(serialize(peptide, notation, **policy), text="A-A")

    monkeypatch.setattr("pyPept.peptide.serialize", incorrect)
    response = client.post(
        "/convert_notation", json={"cabiln": "A-G", "target": "bracket"}
    )
    assert response.status_code == 400
    assert "change the molecular structure" in response.json()["error"]


@pytest.mark.parametrize(
    "route,key", [("/smiles_to_cabiln", "smiles"), ("/to_cabiln", "input")]
)
@pytest.mark.parametrize(
    "smiles,inferred,synthetic",
    [
        ("NCC(=O)O", False, False),
        ("NC(C)C(=O)O", True, False),
        ("N[C@@H](CCC(F)(F)F)C(=O)O", False, True),
    ],
)
def test_conversion_metadata_describes_the_result(
    client, route, key, smiles, inferred, synthetic
):
    response = client.post(route, json={key: smiles})
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["inferred_stereo"] is inferred
    assert bool(data["synthetic_components"]) is synthetic
    assert bool(data["warnings"]) == (inferred or synthetic)


def test_sequence_diagnostics_are_local_without_global_warning_capture(monkeypatch):
    def unexpected_warning(*args, **kwargs):
        raise AssertionError("Request diagnostics leaked into global warnings")

    monkeypatch.setattr("warnings.warn", unexpected_warning)
    messages = []
    Sequence("C.trt(4,2)", warning_sink=messages.append)
    assert any("thioether" in message for message in messages)
    report = Sequence.validate("E.[G(4,1).A(2,1).G(2,1)]")
    assert report.ok
    assert any("bracket" in message for message in report.warnings)


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_comparison_cannot_discard_relative_or_mixture_stereo_groups(group):
    from pyPept.structure import compare_structures

    molecule = "C[C@H](O)[C@H](C)F"
    with pytest.raises(ValueError, match="AND/OR"):
        compare_structures(
            Chem.MolFromSmiles(f"{molecule} |{group}:1,3|"),
            Chem.MolFromSmiles(molecule),
        )


@pytest.mark.parametrize(
    "notation,index,expected",
    [
        ("A-G", 0, "A-K-G"),
        ("!1-A-G-A-!1", 0, "!1-A-K-G-A-!1"),
        ("!1-A-G-A-!1", 2, "!1-A-G-A-K-!1"),
        ("K.{ac(4,2)}-G", 0, "K.{ac(4,2)}-K-G"),
        ("A-G%K-A", 0, "A-K-G%K-A"),
    ],
)
def test_insert_backbone_uses_residue_index(client, notation, index, expected):
    response = client.post(
        "/insert_backbone",
        json={
            "cabiln": notation,
            "after_idx": index,
            "new_abbr": "K",
        },
    )
    assert response.status_code == 200, response.text
    assert response.json()["result"] == expected


@pytest.mark.parametrize(
    "notation,index,host_slot,new_slot,new_abbr,expected",
    [
        ("!1-K-G-A-!1", 0, 4, 2, "ac", "!1-K.{ac(4,2)}-G-A-!1"),
        ("K.{G(4,2)}-A", 2, 1, 2, "ac", "K.{G(4,2).ac(1,2)}-A"),
        ("A-G%K-A", 2, 4, 2, "ac", "A-G%K.{ac(4,2)}-A"),
        ("A-G", 1, 2, 1, "K", "A-G-K"),
    ],
)
def test_insert_bond_targets_selected_residue(
    client, notation, index, host_slot, new_slot, new_abbr, expected
):
    response = client.post(
        "/insert_bond",
        json={
            "cabiln": notation,
            "host_residue_idx": index,
            "new_abbr": new_abbr,
            "r_host": host_slot,
            "r_new": new_slot,
        },
    )
    assert response.status_code == 200, response.text
    assert response.json()["result"] == expected


@pytest.mark.parametrize("index,slot", [(8, 4), (0, 2), (0, 64)])
def test_bond_to_missing_or_occupied_attachment_is_rejected(client, index, slot):
    response = client.post(
        "/insert_bond",
        json={
            "cabiln": "K-A",
            "host_residue_idx": index,
            "new_abbr": "ac",
            "r_host": slot,
            "r_new": 2,
        },
    )
    assert response.status_code == 400
    assert response.json()["error"]


def test_verify_keeps_defined_double_bond_stereo(client):
    expected = canonical("DeltaAbu")
    mol = Chem.MolFromSmiles(expected)
    bond = next(
        b
        for b in mol.GetBonds()
        if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
    )
    bond.SetStereo(
        Chem.BondStereo.STEREOE
        if bond.GetStereo() == Chem.BondStereo.STEREOZ
        else Chem.BondStereo.STEREOZ
    )
    opposite = Chem.MolToSmiles(mol)
    assert opposite != expected
    response = client.post("/verify", json={"cabiln": "DeltaAbu", "smiles": opposite})
    assert response.status_code == 200, response.text
    assert response.json()["match"] is False


@pytest.mark.parametrize(
    "notation,smiles,exact,compatible",
    [
        ("G", "NCC(=O)O", True, True),
        ("A", "NC(C)C(=O)O", False, True),
        ("A", "N[C@H](C)C(=O)O", False, False),
        ("A", "N[13C@@H](C)C(=O)O", False, False),
        ("A", "[NH3+][C@@H](C)C(=O)O", False, False),
        ("G", "NCC(=O)O.[Na+]", False, False),
        ("H", "N[C@@H](Cc1c[nH]cn1)C(=O)O", False, False),
        ("I", "NC([C@@H](C)CC)C(=O)O", False, True),
        ("I", "NC([C@H](C)CC)C(=O)O", False, False),
    ],
)
def test_verify_distinguishes_identity_from_stereo_inference(
    client, notation, smiles, exact, compatible
):
    response = client.post("/verify", json={"cabiln": notation, "smiles": smiles})
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["match"] is exact
    assert data["stereo_compatible"] is compatible
    assert ("warning" in data) == (compatible and not exact)


def test_health_responds_while_rendering(client, monkeypatch):
    started, release = Event(), Event()
    original = rendering._draw_mol

    def slow_draw(*args, **kwargs):
        started.set()
        assert release.wait(10)
        return original(*args, **kwargs)

    monkeypatch.setattr(rendering, "_draw_mol", slow_draw)
    with ThreadPoolExecutor(max_workers=2) as pool:
        render = pool.submit(
            client.post, "/render", json={"cabiln": "G-A-G-A-G", "seed": 4}
        )
        assert started.wait(10)
        try:
            health = pool.submit(client.get, "/health")
            assert health.result(timeout=5).status_code == 200
        finally:
            release.set()
        assert render.result(timeout=10).status_code == 200


@pytest.fixture
def writable_library(tmp_path, monkeypatch):
    path = tmp_path / "monomers.sdf"
    _, existing = monomer_store._load_sdf()
    path.write_text(Chem.SDWriter.GetText(existing["G"], kekulize=False))
    monkeypatch.setattr(monomer_store, "_SDF_PATH", path)
    monomer_store._invalidate_sdf()
    with TestClient(create_app(allow_registration=True)) as client:
        yield path, client
    monomer_store._invalidate_sdf()


def registration_payload(client):
    response = client.post("/preview_monomer", json={"smiles": "N[C@@H](C)C(=O)O"})
    assert response.status_code == 200, response.text
    data = response.json()
    return {key: data[key] for key in ("chuckles", "chem_types", "leaving")} | {
        "abbr": "TestA",
        "name": "Test alanine",
    }


def test_public_app_cannot_modify_library(client):
    assert client.get("/capabilities").json() == {"registration": False}
    response = client.post("/register_monomer", json=registration_payload(client))
    assert response.status_code == 403
    assert client.get("/register").status_code == 403


def test_registration_is_atomic_and_duplicate_safe(writable_library):
    path, client = writable_library
    payload = registration_payload(client)
    with ThreadPoolExecutor(max_workers=2) as pool:
        responses = list(
            pool.map(lambda _: client.post("/register_monomer", json=payload), range(2))
        )
    assert sorted(r.status_code for r in responses) == [200, 400]
    mols = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert len(mols) == 2 and all(mol is not None for mol in mols)
    assert [mol.GetProp("m_abbr") for mol in mols] == ["G", "TestA"]
    assert client.get("/monomer_svg?abbr=TestA").status_code == 200


def test_failed_atomic_replace_keeps_original_file(writable_library, monkeypatch):
    path, client = writable_library
    original = path.read_bytes()
    payload = registration_payload(client)

    def fail_replace(*args):
        raise OSError("Simulated interrupted write")

    monkeypatch.setattr("os.replace", fail_replace)
    response = client.post("/register_monomer", json=payload)
    assert response.status_code == 500
    assert response.json()["request_id"]
    assert "Simulated" not in response.text
    assert path.read_bytes() == original
    assert not list(path.parent.glob(".monomers-*.sdf"))


def test_inconsistent_monomer_metadata_is_rejected(writable_library):
    path, client = writable_library
    original = path.read_bytes()
    payload = registration_payload(client)
    payload["leaving"].pop("1")
    response = client.post("/register_monomer", json=payload)
    assert response.status_code == 400
    assert path.read_bytes() == original


@pytest.mark.parametrize("leaving", ["*", "OC", "[F-]"])
def test_unusable_restored_monomer_is_rejected_before_registration(
    writable_library, leaving
):
    path, client = writable_library
    original = path.read_bytes()
    payload = registration_payload(client)
    payload["leaving"]["2"] = leaving
    response = client.post("/register_monomer", json=payload)
    assert response.status_code == 400, response.text
    assert path.read_bytes() == original
