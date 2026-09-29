"""Canonical formatting and recursive branches through editing and HTTP."""

from dataclasses import replace
from shutil import copyfile

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept.editor import PeptideDocument
from pyPept.inputs import format_source
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.web.app import create_app
from pyPept.web.projects import project_context


@pytest.fixture
def client():
    with TestClient(create_app(allow_registration=False)) as client:
        yield client


def product(source):
    return Chem.MolToSmiles(Molecule(Sequence(source)).get_molecule("ROMol"))


@pytest.mark.parametrize(
    "source",
    [
        "K.[G(4,2).ac(1,2)]-A",
        "K.[G(4,2).[ac(1,2)]]-A",
        "ac-G.!bridge(2,4)%K.!bridge(4,2)-A",
    ],
)
@pytest.mark.parametrize(
    "target,expected",
    [
        ("bracket", "K.[G(4,2).[ac(1,2)]]-A"),
        ("branch", "K.!1(4,2)-A%ac-G.!1(2,4)"),
    ],
)
def test_canonical_route_converges_then_renders_exact_structure(
    client, source, target, expected
):
    response = client.post(
        "/convert_notation",
        json={"cabiln": source, "target": target, "canonical": True},
    )
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["result"] == expected
    assert data["canonical"]["format"] == "cabiln-graph-v1"
    assert data["canonical"]["rdkit"]
    assert set(data["canonical"]["binding"]) == {
        "monomers", "aliases", "reactions", "caps"
    }
    rendered = client.post("/render", json={"cabiln": data["result"]})
    assert rendered.status_code == 200, rendered.text
    molecule = Chem.MolFromMolBlock(rendered.json()["mol_block"])
    assert Chem.MolToSmiles(molecule) == product(source)


@pytest.mark.parametrize("target", ["bracket", "branch"])
def test_existing_layout_conversion_preserves_aliases_and_protection(client, target):
    source = "K.{G(4,2).ac(1,2)}-Ala"
    response = client.post(
        "/convert_notation", json={"cabiln": source, "target": target}
    )
    assert response.status_code == 200, response.text
    assert response.json() == {"result": source, "context": project_context()}


def test_recursive_renderer_keeps_sibling_groups_and_nested_atom_ownership(client):
    source = "K.[K(4,2).[G(1,2).[ac(1,2)]].[ac(4,2)]]-A"
    response = client.post("/render", json={"cabiln": source})
    assert response.status_code == 200, response.text
    data = response.json()
    layout = data["layout"]
    assert layout["segments"][0]["roots"] == [0, 1]
    assert [
        (group["host"], group["roots"], group["parent"], group["members"])
        for group in layout["groups"]
    ] == [
        (0, [2], None, [2, 3, 4, 5]),
        (2, [3], 0, [3, 4]),
        (3, [4], 1, [4]),
        (2, [5], 0, [5]),
    ]
    displayed = layout["segments"][0]["roots"] + [
        index for group in layout["groups"] for index in group["roots"]
    ]
    assert sorted(displayed) == list(range(6))
    ownership = [set(atoms) for atoms in data["residue_map"].values()]
    assert all(ownership)
    assert sum(map(len, ownership)) == len(set.union(*ownership))
    assert len(set.union(*ownership)) == Chem.MolFromMolBlock(
        data["mol_block"]
    ).GetNumAtoms()


def test_editing_nested_repeat_leaves_its_sibling_unchanged():
    source = "K.[K(4,2).[G(1,2)]].[G(1,2)]-A"
    document = PeptideDocument(source)
    first_glycine = next(
        index
        for index, occurrence in enumerate(document.sequence.s_sources)
        if occurrence.token.start == source.index("G(1,2)")
    )
    result = document.attach(document.select(first_glycine), 1, "ac", 2)
    assert result.endswith("].[G(1,2)]-A")
    assert product(result) == product(
        "K.[K(4,2).[G(1,2).[ac(1,2)]]].[G(1,2)]-A"
    )
    assert product(result) != product(
        "K.[K(4,2).[G(1,2)]].[G(1,2).[ac(1,2)]]-A"
    )


@pytest.mark.parametrize(
    "source",
    [
        "K.[K(4,2).[G(1,2)]]-A",
        "K.[K(4,2).G(1,2)]-A",
        "K.[K(4,2)[.G(1,2)]]-A",
    ],
)
def test_new_child_does_not_redirect_existing_children_or_continuations(source):
    document = PeptideDocument(source)
    inner_lysine = next(
        index
        for index, occurrence in enumerate(document.sequence.s_sources)
        if occurrence.token.start == source.index("K(4,2)")
    )
    result = document.attach(document.select(inner_lysine), 4, "ac", 2)
    assert ".{ac(4,2)}" in result
    assert "%" not in result
    assert product(result) == product("K.[K(4,2).[G(1,2)].[ac(4,2)]]-A")


@pytest.mark.parametrize(
    "group",
    [
        "[!r(2,1).G(4,2)]",
        "{!r(2,1).G(4,2)}",
        "[!r(2,1).[G(4,2).G(1,2).ac(1,2)]]",
        "[!r(2,1)[.G(4,2).G(1,2).ac(1,2)]]",
        "[.[!r(2,1)].G(4,2)]",
    ],
)
def test_backbone_insertion_moves_only_marker_from_mixed_recursive_scope(group):
    source = "A.!r(1,2)-K." + group
    document = PeptideDocument(source)
    result = document.insert_backbone(document.select(1), "C")
    arm = (
        "[G(4,2).G(1,2).ac(1,2)]" if "ac(1,2)" in group else "[G(4,2)]"
    )
    assert product(result) == product(f"A.!r(1,2)-K.{arm}-C.!r(2,1)")


@pytest.mark.parametrize("group", ["[.[!r(2,1)]]", "{.[.{!r(2,1)}]}"])
def test_moving_nested_marker_prunes_empty_ancestors(group):
    document = PeptideDocument("A.!r(1,2)-K." + group)
    result = document.insert_backbone(document.select(1), "C")
    assert product(result) == product("A.!r(1,2)-K-C.!r(2,1)")


def test_moving_marker_retains_protection_inherited_from_its_parent():
    source = "K.{.[!r(2,1)].[ac(1,2)].!s(4,2)}%G.!r(1,2)%G.!s(2,4)"
    document = PeptideDocument(source)
    result = document.insert_backbone(document.select(0), "C")
    updated = PeptideDocument(result)
    marker = next(
        marker for marker in updated.peptide.layout.markers
        if marker.label == "!r" and marker.endpoint.slot == 2
    )
    moved = next(
        group for group in updated.peptide.layout.groups if group.id == marker.group
    )
    assert moved.protected
    assert moved.opening == "{"
    assert product(result) == product("ac-K.[G(4,2)]-C-G")


def test_formatting_rejects_a_duplicated_occurrence_mapping(monkeypatch):
    from pyPept import peptide

    serialize = peptide.serialize

    def duplicate(*args, **kwargs):
        return replace(serialize(*args, **kwargs), occurrence_order=(0, 0))

    monkeypatch.setattr(peptide, "serialize", duplicate)
    with pytest.raises(ValueError, match="lose or duplicate a monomer"):
        format_source("G%G", "percent", canonical=True)


def test_export_binding_tracks_content_and_ignores_local_paths(tmp_path, monkeypatch):
    from pyPept import monomer_store

    original = monomer_store.library_path()
    expected = monomer_store.library_binding()
    copied = tmp_path / "monomers.sdf"
    copyfile(original, copied)
    aliases = copied.with_name("monomers.csv")
    copyfile(original.with_name("monomers.csv"), aliases)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(copied))
    assert monomer_store.library_binding() == expected
    with aliases.open("a") as handle:
        handle.write("\n")
    changed = monomer_store.library_binding()
    assert changed["aliases"] != expected["aliases"]
    assert changed["monomers"] == expected["monomers"]
    changed["monomers"] = "caller mutation"
    assert monomer_store.library_binding()["monomers"] == expected["monomers"]
