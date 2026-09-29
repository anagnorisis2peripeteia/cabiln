"""Canonical text must identify the resolved graph, including its decomposition."""

import random
from dataclasses import replace

import pytest
from rdkit import Chem, rdBase

from pyPept.canonical import canonical_convention
from pyPept.molecule import Molecule
from pyPept.peptide import Connection, Endpoint, Peptide, serialize
from pyPept.sequence import Sequence


STYLES = ("percent", "bracket")


def _sequence(source):
    return Sequence(source, track_source=True, warning_sink=lambda _: None)


def _peptide(source):
    return Peptide.from_sequence(_sequence(source))


def _canonical(source, style):
    return serialize(_peptide(source), style, canonical=True).text


def _molecule(source):
    return Chem.MolToSmiles(
        Molecule(_sequence(source), depiction=None).get_molecule(fmt="ROMol")
    )


@pytest.mark.parametrize("style, expected", [
    ("bracket", "K.[G(4,2).[ac(1,2)]]-A"),
    ("percent", "K.!1(4,2)-A%ac-G.!1(2,4)"),
])
@pytest.mark.parametrize("source", [
    "K.[G(4,2).ac(1,2)]-A",
    "K.[G(4,2).[ac(1,2)]]-A",
    "K.[G(4,2)[.ac(1,2)]]-A",
    "ac-G.!bridge(2,4)%Lys.!bridge(4,2)-Ala",
    "K.{G(4,2).ac(1,2)}-A",
])
def test_canonical_example_converges_across_syntax_and_source_layout(
    source, style, expected
):
    assert _canonical(source, style) == expected


@pytest.mark.parametrize("style", STYLES)
def test_aliases_and_sibling_order_converge_but_nested_graphs_remain_distinct(style):
    assert _canonical("Ala-Gly", style) == "A-G"
    siblings = "K.[K(4,2).[K(1,2)].[K(4,2)]]-A"
    swapped = "K.[K(4,2).[K(4,2)].[K(1,2)]]-A"
    deeper = "K.[K(4,2).[K(1,2).[K(4,2)]]]-A"
    assert _canonical(siblings, style) == _canonical(swapped, style)
    assert _canonical(siblings, style) != _canonical(deeper, style)
    changed_slot = "K.[G(1,2).ac(1,2)]-A"
    assert _canonical(changed_slot, style) != _canonical(
        "K.[G(4,2).ac(1,2)]-A", style
    )


GRAPH_CASES = (
    "K.[G(4,2).ac(1,2)]-A",
    "K.[K(4,2).[K(1,2)].[K(4,2)]]-A",
    "K.[K(4,2).[K(1,2).[K(4,2)]]]-A",
    "!ring-G-G-G-G-G-G-!ring",
    "C.!bridge(4,4)-C.!bridge",
    "C.!a(4,4)-C.!b(4,4)-C.!a-C.!b",
    "A-G%A-G%K.[G(4,2).ac(1,2)]-A",
    "ac-C.[TBMB(4,4).!2(5,4).!3(6,4)]-C.!2-C.!3-am",
)


def _shuffle(peptide, seed):
    rng = random.Random(seed)
    nodes = list(peptide.occurrences)
    rng.shuffle(nodes)
    ids = {node.id: 100 + offset * 7 for offset, node in enumerate(nodes)}
    changed = tuple(
        replace(
            node, id=ids[node.id], token="arbitrary_alias", source=None,
            order_key=(rng.randrange(1000000),),
        )
        for node in nodes
    )
    edges = [
        Connection(
            Endpoint(ids[edge.left.occurrence_id], edge.left.slot),
            Endpoint(ids[edge.right.occurrence_id], edge.right.slot),
            f"!source_{seed}_{index}",
        )
        for index, edge in enumerate(peptide.connections)
    ]
    rng.shuffle(edges)
    return Peptide(changed, tuple(edges))


@pytest.mark.parametrize("style", STYLES)
@pytest.mark.parametrize("source", GRAPH_CASES)
def test_graph_permutations_keep_exact_bytes_and_valid_occurrence_maps(source, style):
    peptide = _peptide(source)
    expected = serialize(peptide, style, canonical=True).text
    for seed in range(12):
        permuted = _shuffle(peptide, seed)
        emission = serialize(permuted, style, canonical=True)
        assert emission.text == expected
        assert set(emission.occurrence_order) == {n.id for n in permuted.occurrences}
        assert len(emission.occurrence_order) == len(permuted.occurrences)
        restored = Peptide.from_sequence(
            _sequence(emission.text), emission.occurrence_order
        )
        assert set(restored.connections) == set(permuted.connections)
        assert {n.id: n.definition.key for n in restored.occurrences} == {
            n.id: n.definition.key for n in permuted.occurrences
        }


@pytest.mark.parametrize("source", GRAPH_CASES)
def test_both_writers_are_idempotent_interconvertible_and_preserve_products(source):
    expected = {style: _canonical(source, style) for style in STYLES}
    molecule = _molecule(source)
    for emitted in expected.values():
        assert _molecule(emitted) == molecule
        for style in STYLES:
            assert _canonical(emitted, style) == expected[style]


@pytest.mark.parametrize("style", STYLES)
def test_synthetic_templates_preserve_stereo_isotopes_charge_and_terminal_flags(style):
    template = "[1*]N[C@@H]([13CH2][O-])C([2*:99])=O"
    molecule = Chem.MolFromSmiles(template)
    reordered = Chem.RenumberAtoms(
        molecule, list(reversed(range(molecule.GetNumAtoms())))
    )
    alternate = Chem.MolToSmiles(reordered, canonical=False)
    first = _canonical(f"<{template}>", style)
    assert first == _canonical(f"<{alternate}>", style)
    assert "[13CH2]" in first and "[O-]" in first and "[2*:99]" in first
    assert _molecule(first) == _molecule(f"<{template}>")
    for changed in (
        template.replace("@@", "@"), template.replace("[13CH2]", "C"),
        template.replace("[O-]", "O"), template.replace("[2*:99]", "[2*]"),
    ):
        assert _canonical(f"<{changed}>", style) != first


@pytest.mark.parametrize("style", STYLES)
def test_synthetic_caps_and_recursive_synthetic_branches_keep_their_roles(style):
    source = "<CC(=O)>_-A-_<N>"
    full = "<[2*]C(C)=O>_-A-_<[1*]N>"
    assert _canonical(source, style) == full
    branch = "K.[<[1*]N[C@@H](C[O-])C([2*])=O>(4,2)]-A"
    emitted = _canonical(branch, style)
    assert _molecule(emitted) == _molecule(branch)
    assert _canonical(emitted, style) == emitted
    if style == "bracket":
        assert ".[<" in emitted


@pytest.mark.parametrize("style", STYLES)
def test_canonical_equality_preserves_selected_decomposition_and_named_identity(style):
    combined = "<[1*]NCC(=O)NCC([2*])=O>"
    assert _molecule("G-G") == _molecule(combined)
    assert _canonical("G-G", style) != _canonical(combined, style)
    alanine = _peptide("A").occurrences[0].definition.key[4]
    assert _molecule("A") == _molecule(f"<{alanine}>")
    assert _canonical("A", style) != _canonical(f"<{alanine}>", style)


def test_definition_identity_includes_unused_slot_restoration():
    sequence = _sequence("G")
    before = Peptide.from_sequence(sequence).occurrences[0].definition
    sequence.s_monomers[0]["m_Rgroups"][1] = "[H]"
    changed = Peptide.from_sequence(sequence)
    after = changed.occurrences[0].definition
    assert before.key[4] == after.key[4]
    assert before.key != after.key
    with pytest.raises(ValueError, match="changed a resolved definition"):
        serialize(changed, "percent", canonical=True)


@pytest.mark.parametrize("source, tracking", [("G-G-G", True), ("G-G", False)])
def test_projections_reuse_detached_lazy_definitions_for_equality_and_hash(
    monkeypatch, source, tracking
):
    sequence = Sequence(source, track_source=tracking, warning_sink=lambda _: None)
    expected_template = Chem.MolToSmiles(sequence.s_monomers[0]["m_romol"])
    write_smiles = Chem.MolToSmiles
    calls = []

    def counted(*args, **kwargs):
        calls.append(1)
        return write_smiles(*args, **kwargs)

    monkeypatch.setattr(Chem, "MolToSmiles", counted)
    first = Peptide.from_sequence(sequence)
    second = Peptide.from_sequence(sequence)
    assert not calls
    definition = first.occurrences[0].definition
    assert all(node.definition is definition for node in first.occurrences)
    carbon = next(
        atom for atom in sequence.s_monomers[0]["m_romol"].GetAtoms()
        if atom.GetAtomicNum() == 6
    )
    # Mutate before the first key lookup; a cached key could conceal aliasing.
    carbon.SetIsotope(13)
    assert definition.key[4] == expected_template
    assert first == second
    assert len({first, second}) == 1
    compared_calls = len(calls)
    assert compared_calls > 0
    assert first == second
    assert hash(first) == hash(second)
    assert len(calls) == compared_calls

    changed = Peptide.from_sequence(sequence)
    assert changed.occurrences[0].definition.key != definition.key
    assert changed != first


def test_canonical_layout_has_fresh_spans_and_explicit_parent_ownership():
    source = "K.{K(4,2).[G(1,2).[ac(1,2)]].[ac(4,2)]}-A-A-A"
    peptide = _peptide(source)
    emitted = serialize(peptide, "bracket", canonical=True)
    assert emitted.layout.source == emitted.text
    restored = Peptide.from_sequence(_sequence(emitted.text), emitted.occurrence_order)
    assert emitted.layout == restored.layout
    assert all(group.span is not None for group in emitted.layout.groups)
    assert all(not group.protected for group in emitted.layout.groups)
    assert sum(group.parent is not None for group in emitted.layout.groups) >= 2
    assert serialize(peptide, "preserve").text == source
    assert ".{" in serialize(peptide, "percent").text


def test_canonical_policy_is_explicit_and_versioned():
    peptide = _peptide("Ala-Gly")
    assert serialize(peptide, "percent").text == "Ala-Gly"
    with pytest.raises(ValueError, match="requires 'percent' or 'bracket'"):
        serialize(peptide, "preserve", canonical=True)
    unresolved = Peptide(
        tuple(replace(node, definition=None) for node in peptide.occurrences),
        peptide.connections,
    )
    with pytest.raises(ValueError, match="resolved monomer definitions"):
        serialize(unresolved, "percent", canonical=True)
    assert canonical_convention() == {
        "format": "cabiln-graph-v1",
        "labeling": "rdkit-colored-port-graph-v1",
        "rdkit": rdBase.rdkitVersion,
    }
