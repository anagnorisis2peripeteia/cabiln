"""Recursive scope ownership, complete consumption, and stack-based lowering."""

import pytest
from hypothesis import given, strategies as st
from rdkit import Chem

from _fuzzing import fuzz_settings, record
from _notation_fuzz import branch_cases, edge_set, observed_edges

from pyPept.molecule import Molecule
from pyPept.notation import (
    MAX_BRACKET_DEPTH,
    MAX_NOTATION_CHARACTERS,
    BracketArm,
    bracket_regions,
    normalize_legacy_brackets,
    parse_bracket_group,
    supports_bracket_token,
)
from pyPept.peptide import Peptide, serialize
from pyPept.sequence import Sequence, _expand_inline_caps


def parsed(source):
    return Sequence(source, track_source=True, warning_sink=lambda _: None)


def edges(sequence):
    return {(bond[0], bond[4], bond[2], bond[5]) for bond in sequence.s_bonds}


def product(sequence):
    return Chem.MolToSmiles(Molecule(sequence, depiction=None).mol)


@pytest.mark.parametrize(
    "source",
    [
        "K.[G(4,2).[ac(1,2)]]-A",
        "K.[G(4,2)[.ac(1,2)]]-A",
        "K.{G(4,2).{ac(1,2)}}-A",
        "K.[G(4,2).{ac(1,2)}]-A",
    ],
)
def test_flat_and_explicit_scopes_resolve_the_same_graph(source):
    actual, expected = parsed(source), parsed("K.[G(4,2).ac(1,2)]-A")
    assert (
        edges(actual) == edges(expected) == {(0, 2, 1, 1), (0, 4, 2, 2), (2, 1, 3, 2)}
    )
    assert product(actual) == product(expected)
    assert serialize(Peptide.from_sequence(actual), "preserve").text == source


def test_siblings_and_nested_children_have_different_numbered_connections():
    siblings = parsed("K.[K(4,2).[K(1,2)].[K(4,2)]]-A")
    nested = parsed("K.[K(4,2).[K(1,2).[K(4,2)]]]-A")
    common = {(0, 2, 1, 1), (0, 4, 2, 2), (2, 1, 3, 2)}
    assert edges(siblings) == common | {(2, 4, 4, 2)}
    assert edges(nested) == common | {(3, 4, 4, 2)}
    assert product(siblings) != product(nested)


def test_closing_a_child_restores_the_parent_for_a_flat_step_or_marker():
    flat = parsed("K.[K(4,2).[G(1,2)].ac(4,2)]-A")
    explicit = parsed("K.[K(4,2).[G(1,2)].[ac(4,2)]]-A")
    assert edges(flat) == edges(explicit)
    marker = parsed("K.[K(4,2).[G(1,2)].!x(4,2)]-A%ac.!x(2,4)")
    marker_owner = marker.s_biln.tracker.markers[0].owner
    assert marker_owner == marker.s_sources[3].token
    assert product(marker) == product(flat)


def test_outside_inline_modifiers_remain_independent_of_each_other():
    inline = parsed("K.G(4,2).ac(1,2)-A")
    sequential = parsed("K.[G(4,2).ac(1,2)]-A")
    assert edges(inline) == {(0, 2, 1, 1), (0, 4, 2, 2), (0, 1, 3, 2)}
    assert product(inline) != product(sequential)


@pytest.mark.parametrize(
    "source",
    [
        "K.[[K(4,2).G(1,2)].ac(4,2)]-A",
        "K.[[K(4,2)[.G(1,2)]].ac(4,2)]-A",
        "K.[[[K(4,2).G(1,2)].ac(4,2)]]-A",
    ],
)
def test_legacy_wrappers_keep_their_hub_meaning_without_rewriting_source(source):
    actual = parsed(source)
    expected = parsed("K.[K(4,2).[G(1,2)].[ac(4,2)]]-A")
    assert edges(actual) == edges(expected)
    assert product(actual) == product(expected)
    for monomer, location in zip(actual.s_monomers, actual.s_sources):
        assert source[location.token.start : location.token.end] == monomer["m_abbr"]
    for group in actual.s_biln.tracker.groups:
        assert source[group.span.start : group.span.end]
        if group.legacy:
            assert source[group.span.start] == "."


def test_recursive_legacy_arms_restore_each_parent_cursor():
    source = "K.[K(4,2)[.K(1,2)[.ac(1,2)]][.ac(4,2)]]-A"
    actual = parsed(source)
    expected = parsed("K.[K(4,2).[K(1,2).[ac(1,2)]].[ac(4,2)]]-A")
    assert edges(actual) == edges(expected)
    assert product(actual) == product(expected)


def test_nested_legacy_normalization_uses_the_same_recursive_records():
    source = "K.[G(4,2).[[K(1,2).G(1,2)].ac(4,2)]]-A"
    expected = "K.[G(4,2).[K(1,2)[.G(1,2)][.ac(4,2)]]]-A"
    normalized = normalize_legacy_brackets(source)
    assert normalized == expected
    assert edges(parsed(source)) == edges(parsed(normalized))
    assert normalize_legacy_brackets(normalized) == normalized


@pytest.mark.parametrize("tail", ["", ".ac(4,2)"])
def test_legacy_wrapper_keeps_its_inner_protected_scope_and_exported_hub(tail):
    source = "K.[{K(4,2).G(1,2)}" + tail + "]-A"
    sequence = parsed(source)
    inner = sequence.s_biln.tracker.groups[1]
    assert source[inner.span.start : inner.span.end] == "{K(4,2).G(1,2)}"
    assert inner.protected
    assert [item.protected for item in sequence.s_sources[2:4]] == [True, True]
    if tail:
        assert not sequence.s_sources[4].protected
        assert (2, 4, 4, 2) in edges(sequence)
    normalized = normalize_legacy_brackets(source)
    assert "{K(4,2)[.G(1,2)]}" in normalized
    assert edges(parsed(normalized)) == edges(sequence)
    assert normalize_legacy_brackets(normalized) == normalized


def test_editing_a_legacy_protected_arm_keeps_its_protection():
    from pyPept.editor import PeptideDocument

    document = PeptideDocument("K.[{K(4,2).G(1,2)}]-A")
    result = document.attach(document.select(3), 1, "ac", 2)
    assert "{K(4,2)[.G(1,2).ac(1,2)]}" in result
    assert all(location.protected for location in parsed(result).s_sources[2:])


@pytest.mark.parametrize("partner", [".!x(2,1)", ".!x", ".[.!x]"])
def test_marker_only_scopes_and_bare_partners_keep_their_host(partner):
    actual = parsed("K.[G(4,2).[!x(1,2)]]-A%ac" + partner)
    expected = parsed("K.[G(4,2).ac(1,2)]-A")
    assert product(actual) == product(expected)
    assert len(actual.s_monomers) == 4
    first = actual.s_biln.tracker.markers[0]
    assert first.owner == actual.s_sources[3].token
    assert first.scope == actual.s_biln.tracker.groups[1].span


def test_marker_before_first_monomer_annotates_the_scope_host():
    source = "K.[!x(1,2).[G(4,2)]]-A%ac.!x(2,1)"
    actual = parsed(source)
    assert product(actual) == product(parsed("ac-K.G(4,2)-A"))
    first = actual.s_biln.tracker.markers[0]
    assert first.owner == actual.s_sources[0].token


def test_leading_zero_slot_declarations_compare_by_their_numbered_sites():
    source = "K.[G(04,02).[!x(01,02)]]-A%ac.!x(02,01)"
    assert product(parsed(source)) == product(parsed("K.[G(4,2).ac(1,2)]-A"))


def test_each_recursive_group_records_its_host_parent_and_real_spans():
    source = "K.[K(4,2).{K(1,2).[ac(1,2)]}.[ac(4,2)]]-A"
    sequence = parsed(source)
    locations = sequence.s_sources
    groups = sequence.s_biln.tracker.groups
    assert [group.host for group in groups] == [
        locations[index].token for index in (0, 2, 3, 2)
    ]
    assert [group.parent for group in groups] == [
        None,
        groups[0].span,
        groups[1].span,
        groups[0].span,
    ]
    assert [group.protected for group in groups] == [False, True, True, False]
    assert [location.terminal for location in locations[2:]] == [
        False,
        False,
        True,
        True,
    ]
    assert [location.host for location in locations[2:]] == [
        locations[index].token for index in (0, 2, 3, 2)
    ]
    for location in locations[2:]:
        assert source[location.scope.start : location.scope.end].startswith(".")
    peptide = Peptide.from_sequence(sequence, occurrence_ids=(19, 7, 6, 2, 3, 4))
    assert [group.host for group in peptide.layout.groups] == [19, 6, 2, 6]


def test_implicit_continuation_has_an_owner_without_an_invented_child_group():
    sequence = parsed("K.[G(4,2).ac(1,2)]-A")
    glycine, acetyl = sequence.s_sources[2:]
    assert acetyl.host == glycine.token
    assert acetyl.scope == glycine.scope
    assert acetyl.arm is None
    assert acetyl.parent_scope is None
    assert len(sequence.s_biln.tracker.groups) == 1


def test_recursive_synthetic_templates_are_opaque_and_keep_original_token_spans():
    template = "<[1*]N[C@@H](C[O-])C([2*])=O.[Na+]>"
    ncap = "<CC(=O)[2*]>_"
    ccap = "_<[1*]N>"
    source = f"K.[{template}(4,2).[{ncap}(1,2)]]-A.[{ccap}(2,1)]"
    sequence = parsed(source)
    tokens = [
        source[location.token.start : location.token.end]
        for location in sequence.s_sources
    ]
    assert tokens == ["K", "A", template, ncap, ccap]
    assert all(supports_bracket_token(token) for token in tokens)
    assert serialize(Peptide.from_sequence(sequence), "preserve").text == source
    assert product(sequence)
    raw = ".[" + template + "(4,2).[" + ncap + "(1,2)]]"
    assert list(bracket_regions("K" + raw)) == [(1, len(raw) + 1)]
    assert isinstance(parse_bracket_group(raw).entries[1], BracketArm)
    assert parse_bracket_group(f".[{ccap}(2,1)]").entries[0].token == ccap


def test_parenthesized_library_symbol_can_be_a_recursive_entry_and_group_host():
    symbol = "L_hArg(Et,Et)"
    source = f"{symbol}.[G(1,2).[{symbol}(1,2)]]-A"
    sequence = parsed(source)
    assert [monomer["m_abbr"] for monomer in sequence.s_monomers] == [
        symbol,
        "A",
        "G",
        symbol,
    ]
    assert sequence.s_biln.tracker.groups[0].host == sequence.s_sources[0].token
    assert supports_bracket_token(symbol)
    assert product(sequence)


@pytest.mark.parametrize(
    "source,diagnostic",
    [
        ("K.[G(4,2).[ac(1,2)]-A", "not closed"),
        ("K.[G(4,2).{ac(1,2)]]-A", "mismatched"),
        ("K.[G(4,2).[ac(1,2)garbage]]-A", "unrecognised"),
        ("K.[G(4,2).[ac(1,2)].garbage]-A", "unrecognised"),
        ("K.[G(4,2)..[ac(1,2)]]-A", "unrecognised"),
        ("K.[G(4,2).[]]-A", "no valid"),
        ("K.[G(4,2).[ac(0,2)]]-A", "positive"),
        ("K.[G(4,2).[ac(9,2)]]-A", "R9|rgroup 9"),
        ("K.[G(4,2).[ac(2,2)]]-A", "already|occupied|more than once"),
        ("K.[G(4,2).[!x(1,2)]]-A%ac.!x(2,4)", "R-group conflict"),
        ("K.[G(4,2).[!x(1,2)]]-A", "exactly 2"),
        ("K.[G(4,2).[!x(1,2).!x]]-A%ac.!x", "3rd endpoint"),
    ],
)
def test_every_depth_consumes_and_validates_its_entire_input(source, diagnostic):
    with pytest.raises(ValueError, match=diagnostic):
        parsed(source)
    assert not Sequence.validate(source).ok


def test_lowering_exceeds_python_recursion_depth_without_using_its_call_stack():
    depth = 1200
    source = "K.[G(4,2)" + ".[G(1,2)" * (depth - 1) + "]" * depth + "-A"
    expanded, markers = _expand_inline_caps(source)
    assert expanded.count("G(") == depth
    assert expanded.count(".") == depth
    assert markers == {}


def test_resource_limits_are_explicit_and_actionable():
    source = "K" + ".[G(1,2)" * (MAX_BRACKET_DEPTH + 1)
    source += "]" * (MAX_BRACKET_DEPTH + 1)
    with pytest.raises(ValueError, match="nesting exceeds.*percent segments"):
        _expand_inline_caps(source)
    with pytest.raises(ValueError, match="characters.*smaller documents"):
        Sequence("G" * (MAX_NOTATION_CHARACTERS + 1), track_source=True)


@pytest.mark.fuzz
@fuzz_settings()
@given(case=branch_cases())
def test_fuzz_recursive_spellings_preserve_generated_hosts_and_protection(case):
    expected_edges = edge_set(case.graph.edges)
    record("notation_recursive_input", tokens=case.tokens,
           children=case.children, protected=case.protected)
    canonical = serialize(case.graph.peptide(range(len(case.tokens))),
                          "bracket", canonical=True).text
    for mode in ("explicit", "implied", "legacy", "hub"):
        source, order, protected = case.spelling(mode)
        record("notation_scope:" + mode, nodes=len(order),
               protected=len(protected), source=source)
        sequence = parsed(source)
        peptide = Peptide.from_sequence(sequence, order)
        assert observed_edges(peptide) == expected_edges
        assert {i for i, location in zip(order, sequence.s_sources)
                if location.protected} == protected
        assert serialize(peptide, "preserve").text == source
        normalized = normalize_legacy_brackets(source)
        assert normalize_legacy_brackets(normalized) == normalized
        normalized_graph = Peptide.from_sequence(parsed(normalized), order)
        assert observed_edges(normalized_graph) == expected_edges
        assert serialize(peptide, "bracket", canonical=True).text == canonical
        formatted = serialize(peptide, "percent")
        replay = parsed(formatted.text)
        replay_graph = Peptide.from_sequence(replay, formatted.occurrence_order)
        assert observed_edges(replay_graph) == expected_edges
        assert {i for i, location in zip(formatted.occurrence_order, replay.s_sources)
                if location.protected} == protected


@pytest.mark.fuzz
@fuzz_settings()
@given(depth=st.integers(0, 8), tail=st.sampled_from(("A", "Ala", "Gly")))
def test_fuzz_sibling_and_deeper_children_remain_different(depth, tail):
    left = "K(1,2)" + ".[G(1,2)" * depth + "]" * depth
    siblings = f"K.[K(4,2).[{left}].[K(4,2)]]-{tail}"
    deeper = f"K.[K(4,2).[{left}.[K(4,2)]]]-{tail}"
    common = [(0, 2, 1, 1), (0, 4, 2, 2), (2, 1, 3, 2)]
    common += [(3 + i, 1, 4 + i, 2) for i in range(depth)]
    first = Peptide.from_sequence(parsed(siblings))
    second = Peptide.from_sequence(parsed(deeper))
    record("notation_distinct_scope", depth=depth, siblings=siblings, deeper=deeper)
    assert observed_edges(first) == edge_set(common + [(2, 4, 4 + depth, 2)])
    assert observed_edges(second) == edge_set(common + [(3, 4, 4 + depth, 2)])
    for style in ("percent", "bracket"):
        assert serialize(first, style, canonical=True).text != serialize(
            second, style, canonical=True
        ).text
