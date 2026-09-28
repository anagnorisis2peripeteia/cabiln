"""Compatibility remains explicit while chemistry tests use package imports."""

import importlib.util
import re
from pathlib import Path

import pytest

from pyPept.smiles import smiles_to_cabiln_core


@pytest.fixture
def launcher():
    path = Path(__file__).resolve().parents[1] / "tools" / "live_renderer.py"
    spec = importlib.util.spec_from_file_location("legacy_renderer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_legacy_conversion_import_keeps_the_core_function(launcher):
    assert launcher.smiles_to_cabiln_core is smiles_to_cabiln_core
    source, details = launcher.smiles_to_cabiln_core("NCC(=O)O")
    assert source == "G"
    assert details[0][0] == "G"


def test_legacy_text_helper_retains_unresolved_input_support(launcher):
    source = "ac-C.!1(4,4)-A-A-C.!3(4,4)-am"
    assert launcher._renumber_xlinks(source) == "ac-C.!1(4,4)-A-A-C.!2(4,4)-am"
    assert launcher._renumber_xlinks("X.!27-X.!named-X.!27") == "X.!1-X.!named-X.!1"


def test_legacy_bracket_helper_renumbers_remaining_scaffold_markers(launcher):
    from pyPept.sequence import cabiln_to_bracket

    source = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
    bracket = launcher._renumber_xlinks(cabiln_to_bracket(source))
    assert "TBMB" in bracket
    identifiers = sorted({int(value) for value in re.findall(r"\.!(\d+)", bracket)})
    assert identifiers == list(range(1, len(identifiers) + 1))


def test_legacy_missing_import_remains_attribute_error(launcher):
    with pytest.raises(AttributeError, match="has no attribute"):
        launcher.no_such_legacy_export
