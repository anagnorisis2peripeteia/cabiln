"""Legacy text transforms retain their library-free contract during consolidation."""

import pytest

from pyPept.sequence import biln_to_cabiln, cabiln_to_bracket, cabiln_to_branch


@pytest.mark.parametrize(
    "source,expected",
    [
        ("C(1,3)-G-C(1,3)", "C.!1(4,4)-G-C.!1"),
        ("C(1,3)-G-C(1,03)", "C.!1(4,03)-G-C.!1"),
        ("C(1,03)-G-C(1,03)", "C.!1(03,03)-G-C.!1"),
    ],
)
def test_legacy_slot_mapping_retains_text_spelling_policy(source, expected):
    assert biln_to_cabiln(source) == expected


def test_text_slot_mapping_does_not_require_integer_conversion():
    spelling = "0" * 5000 + "3"
    assert biln_to_cabiln(f"C(1,{spelling})-C(1,3)") == f"C.!1({spelling},4)-C.!1"


@pytest.mark.parametrize(
    "function,source,expected",
    [
        (
            cabiln_to_bracket,
            "Host.!7(4,4)-Tail%Before-Unknown.!7-After",
            "Host.[Unknown(4,4)[.After(2,1)].Before(1,2)]-Tail",
        ),
        (
            cabiln_to_bracket,
            "Host.(4,1)-Tail%Unknown-Other(2,1)",
            "Host.[Unknown(4,1).Other(2,1)]-Tail",
        ),
        (
            cabiln_to_branch,
            "Host.[Unknown(4,2).Other(1,2)]-Tail",
            "Host.!1(4,2)-Tail%Other-Unknown-!1",
        ),
    ],
)
def test_unregistered_symbols_keep_exact_legacy_emission(
    monkeypatch, function, source, expected
):
    import pyPept.sequence as parser

    def no_library(*args, **kwargs):
        raise AssertionError("Text formatting must not resolve a monomer library")

    monkeypatch.setattr(parser, "get_monomer_info", no_library)
    monkeypatch.setattr(parser, "Sequence", no_library)
    assert function(source) == expected


@pytest.mark.parametrize("function", [cabiln_to_bracket, cabiln_to_branch])
@pytest.mark.parametrize(
    "source",
    [
        "Host.{Unknown(4,2).Other(1,2)}-Tail",
        "Host.[Unknown(4,2).Other(1,4)garbage]-Tail",
        "Host.[Unknown(4,2)",
    ],
)
def test_protected_and_malformed_text_keeps_legacy_pass_through(function, source):
    assert function(source) == source
