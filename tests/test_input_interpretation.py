"""Input policies stay consistent across display, CLI and HTTP entry points."""

import pytest
from hypothesis import given, strategies as st
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept.converter import Converter
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence
from pyPept.show import show
from pyPept.web.app import create_app

from _chemistry_fuzz import alternate_smiles
from _fuzzing import fuzz_settings, record


def canonical(molecule):
    return Chem.MolToSmiles(molecule)


@pytest.mark.fuzz
@fuzz_settings()
@given(
    residues=st.lists(st.sampled_from("AGCSV"), min_size=2, max_size=5),
    order=st.integers(min_value=0, max_value=65535),
)
def test_fuzz_supported_imports_agree_with_independent_linear_product(residues, order):
    from pyPept.inputs import detect_input, read_input

    # Literal amino-acid units, independent of library templates and reactions.
    units = {"A": "N[C@@H](C)C(=O)", "G": "NCC(=O)",
             "C": "N[C@@H](CS)C(=O)", "S": "N[C@@H](CO)C(=O)",
             "V": "N[C@@H](C(C)C)C(=O)"}
    source = "".join(units[residue] for residue in residues) + "O"
    expected = canonical(Chem.MolFromSmiles(source))
    inputs = {
        "cabiln": "-".join(residues), "biln": "-".join(residues),
        "fasta": "".join(residues),
        "helm": "PEPTIDE1{" + ".".join(residues) + "}$$$$V2.0",
        "smiles": alternate_smiles(source, order),
    }
    for kind, text in inputs.items():
        record("input." + kind, text=text, order=order)
        parsed = read_input(text, input_format=kind)
        assert parsed.format == kind.upper()
        assert canonical(parsed.assemble()) == expected
        if kind in {"biln", "helm", "smiles"}:
            # C-C is both a valid BILN dipeptide and a valid SMILES alkane.
            # Auto/reference deliberately chooses SMILES; explicit format wins.
            smiles = Chem.MolFromSmiles(text)
            auto_expected = canonical(smiles) if smiles is not None else expected
            assert canonical(detect_input(text).assemble()) == auto_expected


@pytest.mark.fuzz
@fuzz_settings(examples=20)
@given(token=st.sampled_from(("C", "N")), label=st.integers(min_value=1, max_value=999))
def test_fuzz_format_selection_and_legacy_slot_translation_are_explicit(token, label):
    from pyPept.inputs import detect_input, read_input

    record("input.precedence-and-legacy", token=token, label=label)
    assert canonical(read_input(token, input_format="smiles").assemble()) == token
    peptide = read_input(token, input_format="cabiln").assemble()
    assert canonical(peptide) != token
    assert canonical(detect_input(token, policy="display").assemble()) == canonical(peptide)
    assert canonical(detect_input(token, policy="reference").assemble()) == token
    legacy = read_input(f"C({label},3)-A-C({label},3)", input_format="biln")
    modern = read_input(legacy.source, input_format="cabiln")
    assert legacy.sequence.s_bonds == modern.sequence.s_bonds
    assert [edge[4:] for edge in modern.sequence.s_bonds if edge[4:] == [4, 4]] == [[4, 4]]
    with pytest.raises(ValueError, match="CABILN"):
        Converter(biln=legacy.source)


@pytest.mark.fuzz
@fuzz_settings(examples=20)
@given(polymer=st.sampled_from(("RNA", "CHEM")), number=st.integers(min_value=1, max_value=20))
def test_fuzz_non_peptide_helm_and_extended_fasta_are_rejected(polymer, number):
    from pyPept.inputs import read_input

    record("input.unrepresentable", polymer=polymer, number=number)
    with pytest.raises(ValueError):
        read_input(f"{polymer}{number}{{A.G}}$$$$V2.0", input_format="helm")
    # This entry point accepts plain one-letter sequences, not extended tokens.
    with pytest.raises(ValueError):
        read_input(f"[Custom{number}]", input_format="fasta")
    with pytest.raises(ValueError, match="Invalid SMILES"):
        read_input("A" * number, input_format="smiles")


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


@pytest.mark.parametrize(
    "input_format, expected_format, expected_smiles",
    [
        (None, "SMILES", "C"),
        ("auto", "SMILES", "C"),
        ("smiles", "SMILES", "C"),
        ("biln", "BILN", "N[C@@H](CS)C(=O)O"),
        ("cabiln", "CABILN", "N[C@@H](CS)C(=O)O"),
    ],
)
def test_reference_http_respects_selected_input_format(
    input_format, expected_format, expected_smiles
):
    payload = {"input": "C"}
    if input_format is not None:
        payload["input_format"] = input_format
    with TestClient(create_app()) as client:
        response = client.post("/render_reference", json=payload)
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["format"] == expected_format
    assert data["smiles"] == canonical(Chem.MolFromSmiles(expected_smiles))
    assert data["context"]["library_binding"]


@pytest.mark.parametrize("input_format", ["biln", "cabiln"])
@pytest.mark.parametrize("notation", ["percent", "bracket"])
def test_selected_peptide_c_converts_to_cysteine(input_format, notation):
    with TestClient(create_app()) as client:
        response = client.post(
            "/to_cabiln",
            json={"input": "C", "input_format": input_format, "notation": notation},
        )
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["from"] == input_format.upper()
    assert data["cabiln"] == "C"
    assert data["details"] == []
    assert data["context"]["library_binding"]
    product = Molecule(Sequence(data["cabiln"])).get_molecule(fmt="ROMol")
    assert canonical(product) == canonical(Chem.MolFromSmiles("N[C@@H](CS)C(=O)O"))


@pytest.mark.parametrize("input_format", [None, "auto", "smiles"])
def test_methane_conversion_does_not_claim_a_cysteine_residue(input_format):
    payload = {"input": "C"}
    if input_format is not None:
        payload["input_format"] = input_format
    with TestClient(create_app()) as client:
        response = client.post("/to_cabiln", json=payload)
    assert response.status_code == 400, response.text
    assert response.json()["error"]


@pytest.mark.parametrize("route", ["/render_reference", "/to_cabiln"])
def test_explicit_smiles_does_not_fall_back_to_peptide_notation(route):
    with TestClient(create_app()) as client:
        response = client.post(
            route, json={"input": "A", "input_format": "smiles"}
        )
    assert response.status_code == 400, response.text
    assert "Invalid SMILES" in response.json()["error"]


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


@pytest.mark.parametrize(
    "flag, source, expected, depiction",
    [
        (flag, source, expected, depiction)
        for flag, source, expected in [
            ("--biln", "ac-A-G-am", "ac-A-G-am"),
            ("--fasta", "AG", "A-G"),
            (
                "--helm",
                "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0",
                "C.!1(4,4)-C.!1",
            ),
        ]
        for depiction in ("local", "rdkit")
    ] + [("--biln", "ac-A-G-am", "ac-A-G-am", None)],
)
def test_cli_formats_keep_depiction_and_exports(
    flag, source, expected, depiction, monkeypatch, tmp_path
):
    from pyPept.interfaces.run_pyPept import main

    output = tmp_path / "peptide"
    arguments = [
        "run_pyPept", flag, source, "--noconf", "--sdf2D", "--prefix", str(output)
    ]
    if depiction is not None:
        arguments.extend(["--depiction", depiction])
    monkeypatch.setattr("sys.argv", arguments)
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
