"""Source identity and requested topology, beyond mere product buildability."""

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept.editor import EditError, PeptideDocument
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.web.app import app

CHARGED = "<[1*]N[C@@H](C[O-])C([2*])=O>"
RING = "<[1*]N[C@@H](CC%10CCCCC%10)C([2*])=O>"
TRACES = [
    ("ac-A-G-am", ("backbone", 1, "K"), "ac-A-K-G-am"),
    (
        "ac-K.{G(4,2).A(1,2)}-am",
        ("attach", 4, 1, "ac", 2),
        "ac-K.{G(4,2).A(1,2).ac(1,2)}-am",
    ),
    (
        "ac-K.{G(4,2)[.A(1,2)]}-D-am",
        ("connect", 5, 1, 2, 4),
        "ac-K.{G(4,2)[.A(1,2).!1(1,4)]}-D.!1(4,1)-am",
    ),
    ("!1-A-G-K-!1", ("backbone", 0, "C"), "!1-A-C-G-K-!1"),
    ("ac-C-A-C-am", ("connect", 1, 4, 3, 4), "ac-C.!1(4,4)-A-C.!1(4,4)-am"),
    ("A-G%K-A", ("attach", 2, 4, "ac", 2), "A-G%K.{ac(4,2)}-A"),
    ("ac-K.G(4,2)-am%G", ("attach", 3, 1, "ac", 2), "ac-K.G(4,2)-am%ac-G"),
    (f"ac-{CHARGED}-G-am", ("backbone", 1, "K"), f"ac-{CHARGED}-K-G-am"),
]
DISCRIMINATORS = [
    "ac-K.[G(4,2).A(1,2)]-am",
    "ac-K.{G(4,2)[.A(1,2)]}-am",
    "ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am",
    "Ala-Gly",
    f"ac-{RING}-G-am",
    "  A-G  %  K-A  ",
]


def molecule(source):
    return Molecule(Sequence(source)).get_molecule("ROMol")


def smiles(source):
    return Chem.MolToSmiles(molecule(source))


def apply(document, operation):
    name, index, *args = operation
    selected = document.select(index)
    if name == "backbone":
        return document.insert_backbone(selected, *args)
    if name == "attach":
        return document.attach(selected, *args)
    return document.connect(selected, args[0], document.select(args[1]), args[2])


@pytest.mark.parametrize("source,operation,expected", TRACES)
def test_eight_traced_edits(source, operation, expected):
    result = apply(PeptideDocument(source), operation)
    assert result == expected
    # The expected notation was written independently of the editor's emitter.
    # This comparison detects wrong-residue edits even when both products build.
    assert smiles(result) == smiles(expected)


@pytest.mark.parametrize("source", [row[0] for row in TRACES] + DISCRIMINATORS)
def test_provenance_preserves_existing_parser_and_exact_source(source):
    ordinary = Sequence(source)
    document = PeptideDocument(source)
    traced = document.sequence
    assert document.source.encode() == source.encode()
    assert traced.s_biln == ordinary.s_biln
    assert traced.s_chains == ordinary.s_chains
    assert traced.s_bonds == ordinary.s_bonds
    assert [m["m_abbr"] for m in traced.s_monomers] == [
        m["m_abbr"] for m in ordinary.s_monomers
    ]
    assert len(traced.s_sources) == len(traced.s_monomers)
    assert len({o.token for o in traced.s_sources}) == len(traced.s_monomers)
    assert all(source[o.token.start : o.token.end] for o in traced.s_sources)
    assert Chem.MolToSmiles(Molecule(traced).get_molecule("ROMol")) == smiles(source)


def test_exact_occurrence_identity_disambiguates_inline_and_percent_residues():
    source = "ac-K.G(4,2)-am%G"
    document = PeptideDocument(source)
    separate, attached = document.sequence.s_sources[3:]
    assert (
        separate.token.start,
        separate.token.end,
        separate.kind,
        separate.segment,
    ) == (15, 16, "explicit", 1)
    assert (
        attached.token.start,
        attached.token.end,
        attached.kind,
        attached.segment,
    ) == (5, 6, "inline", 0)
    result = document.attach(document.select(3), 1, "ac", 2)
    # Independent RDKit molecule: N-acetyl lysinamide with a glycyl sidechain,
    # plus a separate N-acetyl glycine. The old wrong edit acetylated the sidechain.
    expected = Chem.MolFromSmiles("CC(=O)N[C@@H](CCCCNC(=O)CN)C(N)=O.CC(=O)NCC(=O)O")
    assert Chem.MolToSmiles(molecule(result)) == Chem.MolToSmiles(expected)
    assert smiles(result) != smiles("ac-K.[G(4,2).ac(1,2)]-am%G")


@pytest.mark.parametrize(
    "bracket", ["[G(4,2).A(1,2)]", "{G(4,2).A(1,2)}", "{G(4,2)[.A(1,2)]}"]
)
def test_unselected_protected_and_nested_syntax_is_untouched(bracket):
    source = f"ac-K.{bracket}-am%Ala-Gly"
    document = PeptideDocument(source)
    result = document.insert_backbone(document.select(3), "K")
    assert result == f"ac-K.{bracket}-am%Ala-K-Gly"
    branch = document.sequence.s_sources[5]
    assert source[branch.bracket.start : branch.bracket.end] == "." + bracket
    assert branch.protected is bracket.startswith("{")


def test_nonterminal_branch_attachment_preserves_all_other_edges():
    source = "ac-K.{K(4,2)[.A(1,2)]}-am"
    document = PeptideDocument(source)
    result = document.attach(document.select(3), 4, "ac", 2)
    assert result == "ac-K.{K(4,2).!1(4,2)[.A(1,2)]}-am%ac.!1(2,4)"
    assert smiles(result) == smiles("ac-K.{K(4,2)[.A(1,2)][.ac(4,2)]}-am")


def test_postcondition_rejects_a_buildable_wrong_residue_edit(monkeypatch):
    document = PeptideDocument("ac-K.G(4,2)-am%G")
    append = document._append_to_occurrence
    # Inject the old failure mode: the emitter selects the other identical G.
    monkeypatch.setattr(
        document, "_append_to_occurrence", lambda index, suffix: append(3, suffix)
    )
    with pytest.raises(EditError, match="connections beyond"):
        document.attach(document.select(4), 1, "ac", 2)


def test_selection_cannot_be_reused_after_source_changes():
    old = PeptideDocument("A-G")
    updated = PeptideDocument(old.insert_backbone(old.select(0), "K"))
    with pytest.raises(EditError, match="older source revision"):
        updated.attach(old.select(1), 2, "A", 1)


@pytest.mark.parametrize("token", [CHARGED, RING])
def test_synthetic_source_and_selected_attachment_inspection(token):
    source = token + "-G"
    document = PeptideDocument(source)
    selected = document.sequence.s_sources[0]
    assert source[selected.token.start : selected.token.end] == token
    with TestClient(app) as client:
        rendered = client.post("/render", json={"cabiln": source})
        assert rendered.status_code == 200, rendered.json()
        abbr = document.sequence.s_monomers[0]["m_abbr"]
        inspected = client.get(
            "/monomer_rgroups",
            params={"abbr": abbr, "residue_idx": 0, "cabiln": source},
        )
        assert inspected.status_code == 200, inspected.json()
        slots = {s["slot"]: s for s in inspected.json()["rgroups"]}
        assert not slots[1]["used"]
        assert slots[1]["chem_type"]
        assert slots[2]["used"]
        edited = client.post(
            "/insert_bond",
            json={
                "cabiln": source,
                "host_residue_idx": 0,
                "new_abbr": "ac",
                "r_host": 1,
                "r_new": 2,
            },
        )
        assert edited.status_code == 200, edited.json()
        assert edited.json()["result"] == "ac-" + source
    assert smiles(edited.json()["result"]) == smiles("ac-" + token + "-G")


def test_selected_alias_uses_instance_instead_of_sdf_name_lookup():
    with TestClient(app) as client:
        response = client.get(
            "/monomer_rgroups",
            params={"abbr": "Ala", "residue_idx": 0, "cabiln": "Ala-Gly"},
        )
    assert response.status_code == 200
    assert response.json()["abbr"] == "A"
    assert next(s for s in response.json()["rgroups"] if s["slot"] == 2)["used"]


def test_invalid_selected_index_is_not_silently_treated_as_unbound_tile():
    with TestClient(app) as client:
        response = client.get(
            "/monomer_rgroups", params={"abbr": "G", "residue_idx": 8, "cabiln": "A-G"}
        )
    assert response.status_code == 400
    assert "does not exist" in response.json()["error"]


def test_existing_inline_synthetic_can_be_extended_without_changing_its_token():
    token = "<[1*]NCC([2*])=O>"
    document = PeptideDocument(f"ac-K.{token}(4,2)-am")
    result = document.attach(document.select(3), 1, "ac", 2)
    assert token in result
    assert smiles(result) == smiles("ac-K.{G(4,2).ac(1,2)}-am")


def test_selected_legacy_arm_uses_existing_normalization_before_extension():
    source = "ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am"
    document = PeptideDocument(source)
    result = document.attach(document.select(4), 1, "ac", 2)
    assert result == "ac-K.[K(4,2)[.A(1,2).ac(1,2)][.ac(4,2)]]-am"
    assert smiles(result) == smiles("ac-K.{K(4,2)[.A(1,2).ac(1,2)][.ac(4,2)]}-am")


@pytest.mark.parametrize(
    "source",
    [
        "A.!ring(1,2)-G-A.!ring(2,1)",
        "A.!ring(1,2)-G-A.!ring",
        "!ring-A-G-A-!ring",
    ],
)
def test_backbone_insertion_moves_the_existing_tail_closure(source):
    document = PeptideDocument(source)
    markers = document.sequence.s_biln.tracker.markers
    assert sorted(marker.slot for marker in markers) == [1, 2]
    assert {marker.owner for marker in markers} == {
        document.sequence.s_sources[0].entry,
        document.sequence.s_sources[2].entry,
    }
    result = document.insert_backbone(document.select(2), "K")
    assert smiles(result) == smiles("!ring-A-G-A-K-!ring")


@pytest.mark.parametrize("symbol", ["7New", "_New"])
def test_current_library_names_do_not_need_to_fit_bracket_grammar(
    tmp_path, monkeypatch, symbol
):
    from pyPept import monomer_store
    from pyPept.interfaces.cli_monomer import register_monomer

    _, monomers = monomer_store._load_sdf()
    path = tmp_path / "new-library.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for abbr in ("ac", "am", "K", "G"):
            writer.write(monomers[abbr])
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    try:
        document = PeptideDocument("ac-K-am")
        register_monomer("NCC(=O)O", symbol)
        result = document.attach(document.select(1), 4, symbol, 2)
        assert f"%{symbol}." in result
        assert smiles(result) == smiles("ac-K.G(4,2)-am")
        if symbol.startswith("_"):
            inline = PeptideDocument(f"ac-K.{symbol}(4,2)-am")
            extended = inline.attach(inline.select(3), 1, "ac", 2)
            assert symbol in extended
            assert smiles(extended) == smiles("ac-K.{G(4,2).ac(1,2)}-am")
    finally:
        monomer_store._invalidate_sdf()


@pytest.mark.parametrize("bracket", ["[!r(2,1)]", "{!r(2,1)}"])
def test_backbone_insertion_transfers_a_pure_crosslink_bracket(bracket):
    document = PeptideDocument(f"A.!r(1,2)-G-K.{bracket}")
    result = document.insert_backbone(document.select(2), "C")
    assert result == f"A.!r(1,2)-G-K-C.{bracket}"
    assert smiles(result) == smiles("!r-A-G-K-C-!r")


@pytest.mark.parametrize(
    "bracket,retained,moved",
    [
        ("[!r(2,1).!s(4,4)]", "[!s(4,4)]", "[!r(2,1)]"),
        ("{!s(4,4).!r(2,1)}", "{!s(4,4)}", "{!r(2,1)}"),
    ],
)
def test_transfer_only_outgoing_marker_preserves_other_link_and_bracket_kind(
    bracket, retained, moved
):
    source = f"C.!r(1,2).!s(4,4)-A-C.{bracket}"
    document = PeptideDocument(source)
    result = document.insert_backbone(document.select(2), "G")
    assert result == f"C.!r(1,2).!s(4,4)-A-C.{retained}-G.{moved}"
    assembled = Sequence(result)
    # Independent expected graph: insert G into the cyclic backbone while the
    # two original cysteine R4 sites retain their disulfide connection.
    edges = {frozenset(((b[0], b[4]), (b[2], b[5]))) for b in assembled.s_bonds}
    assert edges == {
        frozenset(((0, 2), (1, 1))),
        frozenset(((1, 2), (2, 1))),
        frozenset(((2, 2), (3, 1))),
        frozenset(((3, 2), (0, 1))),
        frozenset(((0, 4), (2, 4))),
    }
    assert smiles(result) == smiles("!r-C.!s(4,4)-A-C.!s(4,4)-G-!r")


def test_pendant_endpoint_inside_host_token_is_not_the_hosts_endpoint():
    source = "A.!r(1,2)-G-D.[G(4,1).!r(2,1)]-am"
    document = PeptideDocument(source)
    result = document.insert_backbone(document.select(2), "C")
    assert result == "A.!r(1,2)-G-D.[G(4,1).!r(2,1)]-C-am"
    assert smiles(result) == smiles("A.!r(1,2)-G-D.!s(4,1)-C-am%G.!s(1,4).!r(2,1)")
