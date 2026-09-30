"""Bracket lowering must preserve every entry and every declared attachment."""

import pytest
from hypothesis import given, strategies as st
from rdkit import Chem
from unittest.mock import patch

from _fuzzing import fuzz_settings, record

from pyPept.molecule import Molecule
from pyPept.peptide import Peptide, serialize
from pyPept.sequence import (
    Sequence, _expand_inline_caps, cabiln_to_bracket, cabiln_to_branch,
)


def parsed(source):
    return Sequence(source, track_source=True, warning_sink=lambda _: None)


def product(sequence):
    return Chem.MolToSmiles(Molecule(sequence, depiction=None).mol)


@pytest.mark.parametrize(
    "bracket, sequence_message",
    [
        (bracket, None)
        for bracket in [
            "[[K(4,2).A(1,2)garbage].ac(4,2)]",
            "[[?K(4,2).A(1,2)].ac(4,2)]",
            "[[K(4,2).garbage.A(1,2)].ac(4,2)]",
            "[[K(4,2)A(1,2)].ac(4,2)]",
            "[[K(4,2).A(1,2)]garbage.ac(4,2)]",
            "[[K(4,2).A(1,2)].ac(4,2)garbage]",
            "[[K(4,2).A(1,2)]ac(4,2)]",
            "[[K(4,2).A(1,2)]..ac(4,2)]",
            "[[K(4,2)[.A(1,2)]].ac(4,2)garbage]",
            "[[K(4,2)[.A(1,2)garbage]].ac(4,2)]",
            "[[K(4,2).A(1,2)].ac(4,2)",
            "[[K(4,2).A(1,2)}.ac(4,2)]",
            "[K(4,2)[.A(1,2)]garbage.ac(4,2)]",
        ]
    ] + [
        ("[G(4,2)[.bad]]", "sub-bracket|unrecognised"),
        ("[G(4,2)[.A(2,1)garbage]]", "sub-bracket|unrecognised"),
        ("[G(4,2)[.]]", "sub-bracket|unrecognised"),
    ],
)
def test_bracket_parser_rejects_unconsumed_or_unbalanced_regions(bracket, sequence_message):
    source = f"ac-K.{bracket}-am"
    with pytest.raises(ValueError, match="bracket|unrecognised"):
        _expand_inline_caps(source)
    for tracking in (False, True):
        with pytest.raises(ValueError, match=sequence_message):
            Sequence(source, track_source=tracking)
    assert not Sequence.validate(source).ok


@pytest.mark.parametrize("opening,closing", [("[", "]"), ("{", "}")])
@pytest.mark.parametrize("bracket_first", [True, False])
def test_bracket_and_inline_declarations_agree_in_either_order(
    opening, closing, bracket_first
):
    bracket = f"ac-K.{opening}!x(4,1){closing}-am"
    inline = "ac-D.!x(4,4)-am"
    source = "%".join((bracket, inline) if bracket_first else (inline, bracket))
    with pytest.raises(ValueError, match="R-group conflict"):
        parsed(source)


@pytest.mark.parametrize(
    "source",
    [
        "ac-K.[!x(4,1)]-D.!x(4,4)-am",
        "ac-K.!x(4,4)-D.[!x(4,1)]-am",
        "ac-K.[!x(4,1)]-D.[!x(4,4)]-am",
        "ac-K.[G(4,2)[.!x(1,4)]]-D.[!x(4,4)]-am",
        "ac-K.[G(4,2).!x(1,4)]-D.[!x(4,4)]-am",
        "ac-K.[[G(4,2)].!x(1,4)]-D.[!x(4,4)]-am",
    ],
)
def test_all_bracket_positions_share_the_crosslink_declaration_check(source):
    with pytest.raises(ValueError, match="R-group conflict"):
        parsed(source)


@pytest.mark.parametrize(
    "source,explicit",
    [
        ("ac-K.[!x(4,2)]-am%ac-G.!x(2,4)", "ac-K.!x(4,2)-am%ac-G.!x(2,4)"),
        ("ac-K.[!x(4,2)]-am%ac-G.!x", "ac-K.!x(4,2)-am%ac-G.!x(2,4)"),
        ("ac-G.!x%ac-K.[!x(4,2)]-am", "ac-G.!x(2,4)%ac-K.!x(4,2)-am"),
        ("ac-K.[!x(4,2)]-am%ac-G-!x", "ac-K.!x(4,2)-am%ac-G.!x(2,4)"),
        ("ac-G-!x%ac-K.{!x(4,2)}-am", "ac-G.!x(2,4)%ac-K.!x(4,2)-am"),
        ("!x-A-G-!x", "A.!x(1,2)-G.!x(2,1)"),
    ],
)
def test_inverse_bare_and_terminal_partners_keep_their_actual_slots(source, explicit):
    actual, expected = parsed(source), parsed(explicit)
    assert actual.s_bonds == expected.s_bonds
    assert product(actual) == product(expected)


@pytest.mark.parametrize(
    "legacy,modern",
    [
        (
            "ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am",
            "ac-K.[K(4,2)[.A(1,2)][.ac(4,2)]]-am",
        ),
        (
            "ac-K.[[K(4,2)[.A(1,2)]].ac(4,2)]-am",
            "ac-K.[K(4,2)[.A(1,2)][.ac(4,2)]]-am",
        ),
        (
            "ac-K.[[G(4,2)].!x(1,4)]-D.!x(4,1)-am",
            "ac-K.[G(4,2)[.!x(1,4)]]-D.!x(4,1)-am",
        ),
    ],
)
def test_lossless_legacy_lowering_preserves_occurrences_slots_and_source(
    legacy, modern
):
    actual, expected = parsed(legacy), parsed(modern)
    left, right = Peptide.from_sequence(actual), Peptide.from_sequence(expected)
    assert actual.s_bonds == expected.s_bonds
    assert left.connections == right.connections
    assert [(item.id, item.symbol) for item in left.occurrences] == [
        (item.id, item.symbol) for item in right.occurrences
    ]
    assert product(actual) == product(expected)
    assert serialize(left, "preserve").text == legacy
    for occurrence in left.occurrences:
        span = occurrence.source.token
        assert legacy[span.start : span.end] == occurrence.symbol
    assert [(marker.label, marker.endpoint) for marker in left.layout.markers] == [
        (marker.label, marker.endpoint) for marker in right.layout.markers
    ]


def test_synthetic_smiles_brackets_are_protected_before_notation_scanning():
    source = "A-<[1*]N[C@@H](C[O-])C([2*])=O.[Na+]>-G"
    sequence = parsed(source)
    peptide = Peptide.from_sequence(sequence)
    synthetic = peptide.occurrences[1].source.token
    assert source[synthetic.start : synthetic.end] == source[2:-2]
    assert serialize(peptide, "preserve").text == source
    assert product(sequence) == "C[C@H](N)C(=O)N[C@@H](C[O-])C(=O)NCC(=O)O"


@pytest.mark.parametrize(
    "source",
    [
        "ac-K.!1(4,2)-am%ac-G.!1(junk)",
        "ac-K.!1(4,2)-am%ac-G.!1(2,4)(garbage)",
        "ac-K.!1(4,2)-am%ac-G.!1(1,4)",
        "ac-K.!1(4,2)-am%G.!1.!2(1,2)garbage%ac.!2",
        "ac-K.!1(4,2)-am%G.!1.!2(1,2)%ac.!2(1,1)",
        "ac-C.!1(4,4)-C.!2(4,5)-am%TBMB.!1(junk).!2",
        "ac-C.!1(4,4)-C.!2(4,5)-am%TBMB.!1(5,4).!2",
        "ac-K.(4,2)-am%G(junk)-ac(1,2)",
        "ac-K.!1(4,2)-am%ac--G.!1",
        "ac-K.!1(4,2)-am%ac-G.!1.!1",
        "ac-K.!1name(4,2)-am%ac-G.!1name(2,4)",
        "ac-K.!1(4,2)-E.!3(4,4)-am%K.!1.!3.!2(1,2)%ac.!2",
    ],
)
def test_text_branch_conversion_preserves_unsupported_or_conflicting_input(
    monkeypatch, source
):
    import pyPept.sequence as parser

    def no_library(*args, **kwargs):
        raise AssertionError("Text conversion must remain library-free")

    monkeypatch.setattr(parser, "Sequence", no_library)
    monkeypatch.setattr(parser, "get_monomer_info", no_library)
    assert cabiln_to_bracket(source) == source


@pytest.mark.parametrize("partner", [".!2(4,6)", ".[!2(4,6)]", ".{!2(4,6)}"])
def test_text_bracket_conversion_does_not_remove_conflicting_hub_declaration(partner):
    source = f"ac-C.[TBMB(4,4).!2(5,4)]-C{partner}-am"
    with pytest.raises(ValueError, match="R-group conflict"):
        parsed(source)
    assert cabiln_to_branch(source) == source


@pytest.mark.fuzz
@fuzz_settings()
@given(depth=st.integers(1, 8), fault=st.sampled_from((
    "missing_close", "wrong_close", "garbage", "zero_slot", "occupied_slot",
    "missing_slot", "conflicting_pair", "unpaired_tag", "third_endpoint",
)), tag=st.text(alphabet="abcdefghxyz", min_size=1, max_size=8))
def test_fuzz_one_fault_is_rejected_without_graph_repair(depth, fault, tag):
    source = "K.[G(4,2)" + ".[G(1,2)" * (depth - 1) + "]" * depth + "-A"
    if fault in ("conflicting_pair", "unpaired_tag", "third_endpoint"):
        source = f"K.!{tag}(4,2)-A%G.!{tag}(2,4)"
    assert Sequence.validate(source).ok
    if fault == "missing_close":
        position = source.rindex("]")
        broken = source[:position] + source[position + 1:]
    elif fault == "wrong_close":
        position = source.rindex("]")
        broken = source[:position] + "}" + source[position + 1:]
    elif fault == "garbage":
        position = source.index(")") + 1
        broken = source[:position] + "?" + source[position:]
    elif fault == "zero_slot":
        broken = source.replace(",2)", ",0)", 1)
    elif fault == "occupied_slot":
        broken = source.replace("(4,2)", "(2,2)", 1)
    elif fault == "missing_slot":
        broken = source.replace(",2)", ",999)", 1)
    elif fault == "conflicting_pair":
        broken = source.replace(f".!{tag}(2,4)", f".!{tag}(1,4)")
    elif fault == "unpaired_tag":
        broken = source.replace(f"G.!{tag}(2,4)", "G")
    else:
        broken = source + f"%G.!{tag}(2,4)"
    record("notation_malformed:" + fault, original=source, source=broken)
    for tracking in (False, True):
        with pytest.raises(ValueError):
            Sequence(broken, track_source=tracking, warning_sink=lambda _: None)
    assert not Sequence.validate(broken).ok


@pytest.mark.fuzz
@fuzz_settings()
@given(depth=st.integers(1, 8), entries=st.integers(1, 12))
def test_fuzz_parser_limit_boundaries_are_exact(depth, entries):
    import pyPept.notation as grammar

    nested = ".[G(1,2)" * depth + "]" * depth
    flat = ".[" + ".".join(["G(1,2)"] * entries) + "]"
    record("notation_limits", depth=depth, entries=entries, characters=len(flat))
    with patch.object(grammar, "MAX_BRACKET_DEPTH", depth):
        grammar.parse_bracket_group(nested)
        with pytest.raises(ValueError, match="nesting exceeds"):
            grammar.parse_bracket_group(".[G(1,2)" + nested + "]")
    with patch.object(grammar, "MAX_BRACKET_ENTRIES", entries):
        grammar.parse_bracket_group(flat)
        with pytest.raises(ValueError, match="entries"):
            grammar.parse_bracket_group(flat[:-1] + ".G(1,2)]")
    with patch.object(grammar, "MAX_NOTATION_CHARACTERS", len(flat)):
        grammar.parse_bracket_group(flat)
    with patch.object(grammar, "MAX_NOTATION_CHARACTERS", len(flat) - 1):
        with pytest.raises(ValueError, match="characters"):
            grammar.parse_bracket_group(flat)
