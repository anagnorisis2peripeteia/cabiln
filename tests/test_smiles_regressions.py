"""SMILES conversion must preserve molecular structure, not just residue names."""

import pytest
from rdkit import Chem

from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.smiles import (
    StereoInferenceWarning,
    convert_smiles,
    smiles_to_cabiln_core,
)


def _roundtrip(smiles):
    cabiln, details = smiles_to_cabiln_core(smiles)
    rebuilt = Molecule(Sequence(cabiln)).get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(rebuilt, isomericSmiles=True) == Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles), isomericSmiles=True
    )
    return cabiln, details


def test_partly_unspecified_isoleucine_uses_standard_monomer():
    # Retatrutide's source fixture omits Ile alpha stereochemistry. The existing
    # inference policy selects L-Ile; an unused crosslink slot must not select
    # the chemically identical, specialist XlQuat alias instead.
    with pytest.warns(StereoInferenceWarning, match="unspecified"):
        cabiln, _ = smiles_to_cabiln_core("NC([C@@H](C)CC)C(=O)O")
    assert cabiln == "I"


@pytest.mark.parametrize(
    "smiles",
    [
        "N[C@@H](CCC(F)(F)F)C(=O)O",
        "CC(=O)N[C@@H](CCC(F)(F)F)C(=O)N",
        "N[C@H]([C@H](C)CC)C(=O)O",
        "NCC(=O)O.N[C@@H](C)C(=O)O",
        "CC(=O)NCC(=O)N.OC(=O)C",
        "N[C@@H](C/C=C/F)C(=O)O",
        "N[C@@H](C/C=C\\F)C(=O)O",
        "N[13C@@H](C)C(=O)O",
        "CC(=O)N[C@@H](C[NH3+])C(=O)O",
        # Histidine tautomers must not silently become the library's tautomer.
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",
        # Alternative arginine tautomer: matching may normalize it, output may not.
        "N[C@@H](CCCN=C(N)N)C(=O)O",
    ],
)
def test_conversion_preserves_every_atom_and_defined_stereocenter(smiles):
    _roundtrip(smiles)


def test_unsupported_counterion_is_reported_instead_of_discarded():
    with pytest.raises(ValueError, match="component"):
        smiles_to_cabiln_core("NCC(=O)O.[Na+]")


@pytest.mark.parametrize("smiles", ["", " ", "not a SMILES"])
def test_empty_or_invalid_smiles_has_clear_error(smiles):
    with pytest.raises(ValueError, match="SMILES"):
        smiles_to_cabiln_core(smiles)


def test_disconnected_cycles_do_not_share_crosslink_ids():
    parts = [
        Chem.MolToSmiles(Molecule(Sequence(cab)).get_molecule(fmt="ROMol"))
        for cab in ("!1-A-G-!1", "!1-C-A-!1")
    ]
    cabiln, details = _roundtrip(".".join(parts))
    assert "!1" in cabiln and "!2" in cabiln
    assert len(details) == 4


@pytest.mark.parametrize("failure", [KeyboardInterrupt, SystemExit, TypeError])
def test_conversion_does_not_swallow_interrupts_or_programming_errors(
    monkeypatch, failure
):
    def fail(_self, _sequence):
        raise failure("conversion interrupted")

    monkeypatch.setattr(Molecule, "__init__", fail)
    with pytest.raises(failure, match="conversion interrupted"):
        smiles_to_cabiln_core("NCC(=O)O")


@pytest.mark.parametrize(
    "sidechain",
    ["Cc1c[nH]cn1", "CCCN/C(N)=N/[H]"],
    ids=["histidine-tautomer", "arginine-imine-stereo"],
)
def test_local_synthetic_residue_has_the_actual_backbone_attachments(sidechain):
    source = f"N[C@@H]({sidechain})C(=O)O"
    result = convert_smiles(source)
    assert len(result.details) == 1
    assert result.details[0][0] == result.cabiln
    assert result.details[0][1] == 0  # No exact library match is claimed.
    assert result.synthetic_components == (0,)
    assert not result.inferred_stereo
    assert any("local synthetic" in message for message in result.warnings)

    # Both ends must attach at the peptide backbone, rather than another free
    # N on imidazole/guanidine. This tests the converted token in a new context.
    edited = Molecule(Sequence(f"ac-{result.cabiln}-G-am")).get_molecule(fmt="ROMol")
    expected = Chem.MolFromSmiles(f"CC(=O)N[C@@H]({sidechain})C(=O)NCC(N)=O")
    assert Chem.MolToSmiles(edited) == Chem.MolToSmiles(expected)


def test_conversion_status_is_local_to_each_concurrent_call(monkeypatch):
    import warnings
    from concurrent.futures import ThreadPoolExecutor

    def unexpected_warning(*_args, **_kwargs):
        pytest.fail("convert_smiles must return diagnostics, not emit warnings")

    monkeypatch.setattr(warnings, "warn", unexpected_warning)
    inputs = [
        "NCC(=O)O",
        "NC(C)C(=O)O",
        "N[C@@H](C)C(=O)O",
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",
    ]
    with ThreadPoolExecutor(max_workers=4) as pool:
        glycine, inferred, alanine, synthetic = list(pool.map(convert_smiles, inputs))
    assert glycine.cabiln == "G"
    assert glycine.warnings == ()
    assert not glycine.inferred_stereo
    assert inferred.cabiln == "A"
    assert inferred.inferred_stereo
    assert any("unspecified" in message for message in inferred.warnings)
    assert alanine.cabiln == "A"
    assert alanine.warnings == ()
    assert not alanine.inferred_stereo
    assert synthetic.synthetic_components == (0,)
    assert not synthetic.inferred_stereo
    assert all("unspecified" not in message for message in synthetic.warnings)


def test_accepted_assembly_diagnostics_are_returned_without_warning_side_effects(
    monkeypatch,
):
    import warnings

    molecule = Molecule(
        Sequence("ac-C.acm(4,2)-G-am", warning_sink=lambda _message: None)
    ).get_molecule(fmt="ROMol")

    def unexpected_warning(*_args, **_kwargs):
        pytest.fail("assembly diagnostics must use the request-local sink")

    monkeypatch.setattr(warnings, "warn", unexpected_warning)
    result = convert_smiles(Chem.MolToSmiles(molecule))
    assert result.cabiln == "ac-C.acm(4,2)-G-am"
    assert any("thioether" in message for message in result.warnings)
    assert result.synthetic_components == ()


def test_synthetic_component_indices_refer_to_source_components():
    result = convert_smiles("NCC(=O)O.N[C@@H](Cc1c[nH]cn1)C(=O)O")
    assert result.synthetic_components == (1,)
    assert result.cabiln.startswith("G%<")
    assert len(result.details) == 2
    assert result.warnings[0].startswith("Component 2:")


def test_coarse_fallback_does_not_claim_editable_residue_decomposition():
    result = convert_smiles("CC(=O)N[C@@H](CCC(F)(F)F)C(=O)N")
    assert result.details == []
    assert result.synthetic_components == (0,)
    assert any("coarse synthetic" in message for message in result.warnings)
    assert any("residue boundaries" in message for message in result.warnings)


def test_conversion_follows_library_selection_and_new_registrations(
    tmp_path, monkeypatch
):
    import os

    from pyPept.interfaces.cli_monomer import register_monomer

    first, second = tmp_path / "first.sdf", tmp_path / "second.sdf"
    register_monomer("NCC(=O)O", symbol="CustomOne", sdf_path=first)
    second.write_bytes(first.read_bytes().replace(b"CustomOne", b"CustomTwo"))
    stamp = first.stat()
    os.utime(second, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))

    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(first))
    assert convert_smiles("NCC(=O)O").cabiln == "CustomOne"
    # Equal size/timestamp at another path must not reuse the old name.
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(second))
    assert convert_smiles("NCC(=O)O").cabiln == "CustomTwo"
    register_monomer("N[C@@H](C)C(=O)O", symbol="FreshAla")
    source = "NCC(=O)N[C@@H](C)C(=O)O"
    result = convert_smiles(source)
    assert result.cabiln == "CustomTwo-FreshAla"
    assert result.synthetic_components == ()
    rebuilt = Molecule(Sequence(result.cabiln)).get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(Chem.MolFromSmiles(source))
