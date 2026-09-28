"""Input policies stay consistent across display, CLI and HTTP entry points."""

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept.converter import Converter
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.show import show
from pyPept.web.app import create_app


def canonical(molecule):
    return Chem.MolToSmiles(molecule)


@pytest.mark.parametrize(
    "source",
    ["C.!1(4,4)-A-C.!1", "K.[G(4,2).ac(1,2)]-A"],
)
def test_show_exports_valid_cabiln_with_attachment_punctuation(source, tmp_path):
    output = tmp_path / "peptide.mol"
    show(source, fmt="mol", output=output, open_after=False)
    assert output.exists()
    exported = Chem.MolFromMolFile(str(output))
    expected = Molecule(Sequence(source)).get_molecule(fmt="ROMol")
    assert canonical(exported) == canonical(expected)


def test_helm_generated_cabiln_keeps_sidechain_slots_in_http_import():
    helm = "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0"
    emitted = Converter(helm=helm).get_biln()
    assert emitted == "C.!1(4,4)-C.!1"
    # This public input remains deliberately legacy-only, despite its output name.
    with pytest.raises(ValueError, match="CABILN"):
        Converter(biln=emitted)
    expected = canonical(Molecule(Sequence(emitted)).get_molecule(fmt="ROMol"))
    with TestClient(create_app()) as client:
        rendered = client.post("/render_reference", json={"input": helm})
        assert rendered.status_code == 200, rendered.text
        assert rendered.json()["format"] == "HELM"
        assert rendered.json()["smiles"] == expected
        for notation in ("percent", "bracket"):
            converted = client.post(
                "/to_cabiln", json={"input": helm, "notation": notation}
            )
            assert converted.status_code == 200, converted.text
            text = converted.json()["cabiln"]
            assert (
                canonical(Molecule(Sequence(text)).get_molecule(fmt="ROMol"))
                == expected
            )


@pytest.mark.parametrize(
    "argument, expected",
    [
        ({"biln": "C(7,3)-A-C(7,3)"}, "C.!1(4,4)-A-C.!1"),
        ({"chuckles": "C(7,3).A.C(7,3)"}, "C.!1(4,4)-A-C.!1"),
        ({"biln": "Unregistered(2,3)-Other(2,1)"}, "Unregistered.!1(4,1)-Other.!1"),
    ],
)
def test_converter_text_only_inputs_keep_legacy_slot_meaning(argument, expected):
    converted = Converter(**argument)
    assert converted.get_biln() == expected
    assert converted.polymerinfo["bonds"][0][2] == 3
    assert ":R3-" in converted.get_helm()


@pytest.mark.parametrize(
    "source, expected",
    [("C", "C"), ("N", "N"), ("A-G", "A-G")],
)
def test_display_preserves_peptide_precedence_for_bare_tokens(source, expected):
    from pyPept.show import _to_rdmol

    molecule, count = _to_rdmol(source)
    sequence = Sequence(expected)
    assert count == sequence.length()
    assert canonical(molecule) == canonical(
        Molecule(sequence).get_molecule(fmt="ROMol")
    )


@pytest.mark.parametrize("source", ["C", "N", "[C]", "C1CC1", "NCC(=O)O"])
def test_reference_preserves_smiles_precedence(source):
    from pyPept.inputs import detect_input

    parsed = detect_input(source)
    assert parsed.format == "SMILES"
    assert canonical(parsed.molecule) == canonical(Chem.MolFromSmiles(source))


def test_explicit_legacy_input_maps_r3_once_and_keeps_modern_r3_distinct():
    from pyPept.inputs import read_input

    legacy = read_input("C(9,3)-A-C(9,3)", input_format="biln", track_source=True)
    assert legacy.source == "C.!9(4,4)-A-C.!9"
    assert any(edge[4:] == [4, 4] for edge in legacy.sequence.s_bonds)
    modern = read_input("ac-D.!1(4,3)-A.!1(3,4)-G-am")
    expected = Molecule(Sequence("ac-D.!1(4,3)-A.!1(3,4)-G-am")).get_molecule(
        fmt="ROMol"
    )
    assert canonical(modern.assemble()) == canonical(expected)


@pytest.mark.parametrize("depiction", ["local", "rdkit"])
@pytest.mark.parametrize(
    "flag, source, expected",
    [
        ("--biln", "ac-A-G-am", "ac-A-G-am"),
        ("--fasta", "AG", "A-G"),
        (
            "--helm",
            "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0",
            "C.!1(4,4)-C.!1",
        ),
    ],
)
def test_cli_formats_keep_depiction_and_exports(
    flag, source, expected, depiction, monkeypatch, tmp_path
):
    from pyPept.interfaces.run_pyPept import main

    output = tmp_path / "peptide"
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_pyPept",
            flag,
            source,
            "--depiction",
            depiction,
            "--noconf",
            "--sdf2D",
            "--prefix",
            str(output),
        ],
    )
    main()
    assert output.with_suffix(".png").read_bytes().startswith(b"\x89PNG")
    exported = Chem.MolFromMolFile(str(output.with_suffix(".sdf")))
    assert canonical(exported) == canonical(
        Molecule(Sequence(expected)).get_molecule(fmt="ROMol")
    )


def test_format_verification_assembles_each_representation_once(monkeypatch):
    from pyPept.inputs import convert_input

    original = Molecule.__init__
    sources = []

    def assemble(self, sequence, *args, **kwargs):
        sources.append(sequence.s_inputbiln)
        return original(self, sequence, *args, **kwargs)

    monkeypatch.setattr(Molecule, "__init__", assemble)
    kind, result = convert_input("K.[G(4,2).ac(1,2)]-A", "percent")
    assert kind == "BILN"
    assert len(sources) == 2
    assert canonical(Molecule(Sequence(result)).get_molecule(fmt="ROMol")) == canonical(
        Molecule(Sequence("K.[G(4,2).ac(1,2)]-A")).get_molecule(fmt="ROMol")
    )
