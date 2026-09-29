"""Bracket lowering must preserve every entry and every declared attachment."""

import pytest
from rdkit import Chem

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
    "bracket",
    [
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
    ],
)
def test_bracket_parser_rejects_unconsumed_or_unbalanced_regions(bracket):
    source = f"ac-K.{bracket}-am"
    with pytest.raises(ValueError, match="bracket|unrecognised"):
        _expand_inline_caps(source)
    for tracking in (False, True):
        with pytest.raises(ValueError):
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
