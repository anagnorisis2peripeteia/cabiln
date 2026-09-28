"""Behavior regressions from the core correctness audit."""

import pytest
from rdkit import Chem

from pyPept.conformer import Conformer
from pyPept.converter import Converter
from pyPept.interfaces.monomer_pipeline import pre_activate
from pyPept.molecule import Molecule
from pyPept.monomerlib import MonomerLibrary
from pyPept.sequence import Sequence, cabiln_to_bracket, cabiln_to_branch


def smiles(notation):
    return Chem.MolToSmiles(Molecule(Sequence(notation)).mol)


@pytest.mark.parametrize("notation", ["A-ac-G", "am-A", "A-fmoc-G"])
def test_dash_requires_left_r2_and_right_r1(notation):
    with pytest.raises(ValueError, match="backbone|R2"):
        Sequence(notation)
    assert not Sequence.validate(notation).ok


@pytest.mark.parametrize("notation", ["A-<notASmiles>-G", "<broken>_-A", "A-_<broken>"])
def test_invalid_synthetic_smiles_is_not_replaced(notation):
    with pytest.raises(ValueError, match="Synthetic"):
        Sequence(notation)
    assert not Sequence.validate(notation).ok


@pytest.mark.parametrize("notation", [None, 12, "", " ", "A--G", "unknownmonomer"])
def test_invalid_input_raises_value_error_and_reports_failure(notation):
    with pytest.raises(ValueError):
        Sequence(notation)
    assert not Sequence.validate(notation).ok


@pytest.mark.parametrize(
    "notation",
    [
        "ac-K.[G(4,2)[.bad]]-am",
        "ac-K.[G(4,2)[.A(2,1)garbage]]-am",
        "ac-K.[G(4,2)[.]]-am",
    ],
)
def test_sub_brackets_do_not_drop_invalid_content(notation):
    with pytest.raises(ValueError, match="sub-bracket|unrecognised"):
        Sequence(notation)


@pytest.mark.parametrize(
    "notation",
    [
        "ac-C.!1(4,4).!1(4,4)-am",
        "ac-C.!1(4,4).!2(4,4)-am%TBMB.!1.!2",
        "ac-C.!1(4,4).!2(4,4)-C.!1-C.!2-am",
        "ac-G.!1(1,2)-G.!1-am",
    ],
)
def test_each_attachment_slot_is_consumed_once(notation):
    with pytest.raises(ValueError, match="already|same attachment"):
        Sequence(notation)


def test_validate_exposes_backbone_slots():
    assert Sequence.validate("A-G").bonds == [(0, 2, 1, 1)]


@pytest.mark.parametrize("tag", ["bridge", "a2", "1"])
def test_named_crosslinks_assemble(tag):
    assert smiles(f"ac-C.!{tag}(4,4)-G-C.!{tag}-am") == smiles("ac-C.!1(4,4)-G-C.!1-am")


@pytest.mark.parametrize(
    "notation",
    [
        "A-G%K-L",
        "ac-A-am%ac-C.!1(4,4)-C.!1-am",
        "ac-D.!1(4,1)-am%G.!1-am%K-L",
    ],
)
def test_branch_conversion_preserves_unattached_segments(notation):
    assert smiles(cabiln_to_bracket(notation)) == smiles(notation)


def test_midpoint_branch_uses_independent_arms():
    notation = "ac-C.!1(4,4)-am%G-C.!1-A"
    converted = cabiln_to_bracket(notation)
    assert smiles(converted) == smiles(notation)


def test_protected_bracket_survives_conversion():
    notation = "ac-K.{AEEA(4,2).ac(1,2)}-G-am"
    assert cabiln_to_branch(notation) == notation
    assert smiles(notation) == smiles(notation.replace("{", "[").replace("}", "]"))


@pytest.mark.parametrize(
    "notation",
    [
        "K.[E_g(4,4).L(2,1).A(1,2)]-G-am",
        "K.[G(4,2).A(1,2)garbage]-am",
        "K.[!1(4,4).!2(5,4)]-D.!1-E.!2",
    ],
)
def test_branch_conversion_does_not_invent_or_discard_bracket_entries(notation):
    assert cabiln_to_branch(notation) == notation


def test_branch_conversion_preserves_monomer_before_nested_crosslink():
    notation = "ac-K.[E_g(4,4).G(1,2)[.!2(1,4)]]-D.!2-am"
    assert smiles(cabiln_to_branch(notation)) == smiles(notation)


def test_bracket_conversion_preserves_inline_cap_inside_branch():
    notation = "ac-D.!1(4,1)-am%G.!1-K.boc(4,2)-am"
    assert smiles(cabiln_to_bracket(notation)) == smiles(notation)


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("CN[C@@H](CS)C(=O)O", "thiol"),
        ("OC(=O)[C@@H]1C[C@@H](O)CN1", "hydroxyl"),
    ],
)
def test_sidechains_never_use_reserved_r3(raw, expected):
    result = pre_activate(raw)
    assert 3 not in result.chem_types
    assert result.chem_types[4] == expected


def test_library_chuckles_returns_parseable_structure():
    library = MonomerLibrary()
    encoded = library.chuckles_smiles(library.monomer_df.loc["A"])
    assert encoded
    assert Chem.MolFromSmiles(encoded) is not None


@pytest.mark.parametrize(
    "notation, expected",
    [
        ("ac-C.trt(4,2)-G-am", "CG"),
        ("!1-A-G-K-!1", "AGK"),
        ("ac-C.!1(4,4)-G-C.!1-am", "CGC"),
        ("ac-K.[AEEA(4,2).ac(1,2)]-G-am", "KG"),
        ("ac-Aib-G-am", "AG"),
        ("A-<C>-G", "AAG"),
    ],
)
def test_conformer_main_chain_projection_uses_cabiln_parser(notation, expected):
    assert Conformer.get_peptide(notation) == expected


@pytest.mark.parametrize("sidechain", ["C[N+](C)(C)C", "C[O-]", "C%10CCCCC%10"])
def test_synthetic_charge_stereo_and_ring_labels_survive_assembly(sidechain):
    notation = f"<[1*]N[C@@H]({sidechain})C([2*])=O>"
    expected = Chem.MolToSmiles(Chem.MolFromSmiles(f"N[C@@H]({sidechain})C(=O)O"))
    assert smiles(notation) == expected


def test_terminal_restoration_preserves_leaving_atom_charge_and_isotope():
    sequence = Sequence("A")
    sequence.s_monomers[0]["m_Rgroups"][1] = "[18O-]"
    molecule = Molecule(sequence).mol
    expected = Chem.MolFromSmiles("N[C@@H](C)C(=O)[18O-]")
    assert Chem.MolToSmiles(molecule) == Chem.MolToSmiles(expected)


def test_terminal_restoration_rejects_unsupported_multiatom_leaving_group():
    sequence = Sequence("A")
    sequence.s_monomers[0]["m_Rgroups"][1] = "OC"
    with pytest.raises(ValueError, match="leaving group"):
        Molecule(sequence)


@pytest.mark.parametrize("helm", ["", "invalid", "PEPTIDE1{A}$$$$V2.0$extra"])
def test_invalid_helm_is_a_value_error(helm):
    with pytest.raises(ValueError, match="HELM"):
        Converter(helm=helm)


def test_helm_connection_cannot_use_residue_zero():
    with pytest.raises(ValueError, match="residue"):
        Converter(helm="PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,0:R3-2:R3$$$V2.0")


@pytest.mark.parametrize("argument", [{"biln": "C(1,3)-A"}, {"chuckles": "C(1,3).A"}])
def test_incomplete_legacy_bond_is_not_silently_dropped(argument):
    with pytest.raises(ValueError, match="bond"):
        Converter(**argument)


def test_legacy_converter_rejects_cabiln_instead_of_splitting_inline_dots():
    with pytest.raises(ValueError, match="CABILN"):
        Converter(biln="C.!1(4,4)-A-C.!1")


def test_validator_cli_accepts_cabiln(monkeypatch, capsys):
    from pyPept.interfaces.cli_BILN_validityChecker import run_as_standalone

    monkeypatch.setattr("sys.argv", ["validate", "--biln", "ac-C.trt(4,2)-G-am"])
    run_as_standalone()
    assert "NONE\tTrue\tTrue" in capsys.readouterr().out


def test_validator_cli_fails_for_missing_backbone_attachment(monkeypatch, capsys):
    from pyPept.interfaces.cli_BILN_validityChecker import run_as_standalone

    monkeypatch.setattr("sys.argv", ["validate", "--biln", "A-ac-G"])
    with pytest.raises(SystemExit) as error:
        run_as_standalone()
    assert error.value.code == 1
    assert "backbone" in capsys.readouterr().out


def test_structure_cli_can_write_2d_without_pdb_metadata(monkeypatch, tmp_path):
    from pyPept.interfaces.run_pyPept import main

    prefix = tmp_path / "peptide"
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_pyPept",
            "--biln",
            "ac-A-G-am",
            "--noconf",
            "--sdf2D",
            "--prefix",
            str(prefix),
        ],
    )
    main()
    assert prefix.with_suffix(".png").stat().st_size > 0
    output = Chem.MolFromMolFile(str(prefix.with_suffix(".sdf")))
    assert Chem.MolToSmiles(output) == smiles("ac-A-G-am")


def test_pdb_naming_reports_missing_library_metadata():
    from pyPept.sequence import correct_pdb_atoms

    with pytest.raises(ValueError, match="pdbName.*--noconf"):
        correct_pdb_atoms(Sequence("ac-A-G-am"))


@pytest.mark.parametrize("preexisting", [False, True])
def test_registration_supports_new_and_empty_custom_libraries(tmp_path, preexisting):
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    if preexisting:
        path.touch()
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    loaded = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert len(loaded) == 1 and loaded[0].GetProp("symbol") == "CustomGly"


def test_cli_registration_rejects_duplicate_symbol_without_changing_file(tmp_path):
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="CustomGly.*already exists"):
        register_monomer("N[C@@H](C)C(=O)O", symbol="CustomGly", sdf_path=path)
    assert path.read_bytes() == original


@pytest.mark.parametrize("collision_property", ["m_abbr", "symbol"])
@pytest.mark.parametrize("new_property", ["m_abbr", "symbol"])
def test_registration_checks_both_monomer_identifiers(
    tmp_path,
    collision_property,
    new_property,
):
    from pyPept import monomer_store

    path = tmp_path / "custom.sdf"
    _, by_abbr = monomer_store._load_sdf()
    existing = Chem.Mol(by_abbr["G"])
    existing.SetProp("symbol", "ExistingSymbol")
    existing.SetProp("m_abbr", "ExistingAbbreviation")
    path.write_text(Chem.SDWriter.GetText(existing, kekulize=False))
    original = path.read_bytes()
    incoming = Chem.Mol(by_abbr["A"])
    incoming.SetProp("symbol", "NewSymbol")
    incoming.SetProp("m_abbr", "NewAbbreviation")
    incoming.SetProp(new_property, existing.GetProp(collision_property))
    with pytest.raises(ValueError, match="already exists"):
        monomer_store.register_molecule(incoming, sdf_path=path)
    assert path.read_bytes() == original


def test_registration_rejects_active_csv_alias_but_allows_inactive_alias(tmp_path):
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    (tmp_path / "monomers.csv").write_text(
        "token,synonyms\nCustomGly,GlyAlias\nMissingToken,UnusedAlias\n"
    )
    original = path.read_bytes()
    with pytest.raises(ValueError, match="GlyAlias.*already exists"):
        register_monomer("N[C@@H](C)C(=O)O", symbol="GlyAlias", sdf_path=path)
    assert path.read_bytes() == original
    register_monomer("N[C@@H](C)C(=O)O", symbol="UnusedAlias", sdf_path=path)
    assert len(Chem.SDMolSupplier(str(path))) == 2


def test_registration_cannot_shadow_a_paired_cap_alias(tmp_path):
    from pyPept import monomer_store
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    _, by_abbr = monomer_store._load_sdf()
    path.write_text(
        "".join(
            Chem.SDWriter.GetText(by_abbr[symbol], kekulize=False)
            for symbol in ("Bn_", "_Bn")
        )
    )
    original = path.read_bytes()
    with pytest.raises(ValueError, match="Bn.*already exists"):
        register_monomer("NCC(=O)O", symbol="Bn", sdf_path=path)
    assert path.read_bytes() == original


def test_concurrent_cli_registration_has_one_winner(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"

    def register_once(_):
        try:
            register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
        except ValueError as exc:
            assert "already exists" in str(exc)
            return False
        return True

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(register_once, range(2)))
    assert sorted(results) == [False, True]
    loaded = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert len(loaded) == 1 and loaded[0] is not None


@pytest.mark.parametrize("preexisting", [False, True])
def test_failed_cli_registration_keeps_library_unchanged(
    tmp_path,
    monkeypatch,
    preexisting,
):
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    if preexisting:
        register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    original = path.read_bytes() if preexisting else None

    def fail_replace(*args):
        raise OSError("Simulated interrupted write")

    monkeypatch.setattr("os.replace", fail_replace)
    with pytest.raises(OSError, match="Simulated interrupted write"):
        register_monomer("N[C@@H](C)C(=O)O", symbol="CustomAla", sdf_path=path)
    if preexisting:
        assert path.read_bytes() == original
    else:
        assert not path.exists()
    assert not list(tmp_path.glob(".monomers-*.sdf"))


def test_cli_registration_updates_default_library_reads(tmp_path, monkeypatch):
    from pyPept import monomer_store
    from pyPept.interfaces import cli_monomer

    path = tmp_path / "monomers.sdf"
    cli_monomer.register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    monkeypatch.setattr(monomer_store, "_SDF_PATH", path)
    monkeypatch.setattr(cli_monomer, "_default_sdf_path", lambda: path)
    monomer_store._invalidate_sdf()
    try:
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly"}
        cli_monomer.register_monomer("N[C@@H](C)C(=O)O", symbol="CustomAla")
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly", "CustomAla"}
    finally:
        monomer_store._invalidate_sdf()


def test_monomer_cli_duplicate_is_an_actionable_failure(tmp_path, monkeypatch, capsys):
    from pyPept.interfaces import cli_monomer

    path = tmp_path / "custom.sdf"
    cli_monomer.register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    monkeypatch.setattr(
        "sys.argv",
        [
            "pyPept-monomer-add",
            "--smiles",
            "N[C@@H](C)C(=O)O",
            "--symbol",
            "CustomGly",
            "--sdf",
            str(path),
        ],
    )
    with pytest.raises(SystemExit) as error:
        cli_monomer.main()
    assert error.value.code == 1
    assert "CustomGly' already exists" in capsys.readouterr().err


def test_custom_library_is_shared_by_registration_and_core_readers(
    tmp_path, monkeypatch
):
    from pyPept import monomer_store
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    try:
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly"}
        assert set(MonomerLibrary().GetMonomerNames()) == {"CustomGly"}
        register_monomer("N[C@@H](C)C(=O)O", symbol="CustomAla")
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly", "CustomAla"}
        expected = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
        assert smiles("CustomGly-CustomAla") == Chem.MolToSmiles(expected)
    finally:
        monomer_store._invalidate_sdf()


def test_library_cache_distinguishes_paths_with_identical_timestamps(
    tmp_path, monkeypatch
):
    import os

    from pyPept import monomer_store
    from pyPept.interfaces.cli_monomer import register_monomer

    first, second = tmp_path / "first.sdf", tmp_path / "second.sdf"
    register_monomer("NCC(=O)O", symbol="CustomOne", sdf_path=first)
    second.write_bytes(first.read_bytes().replace(b"CustomOne", b"CustomTwo"))
    stat = first.stat()
    os.utime(second, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    try:
        monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(first))
        assert set(monomer_store._load_sdf()[1]) == {"CustomOne"}
        monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(second))
        assert set(monomer_store._load_sdf()[1]) == {"CustomTwo"}
    finally:
        monomer_store._invalidate_sdf()


def test_library_readers_observe_external_file_updates(tmp_path, monkeypatch):
    from pyPept import monomer_store
    from pyPept.interfaces.cli_monomer import register_monomer

    path, extra = tmp_path / "custom.sdf", tmp_path / "extra.sdf"
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=path)
    register_monomer("N[C@@H](C)C(=O)O", symbol="CustomAla", sdf_path=extra)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    try:
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly"}
        path.write_bytes(path.read_bytes() + extra.read_bytes())
        assert set(monomer_store._load_sdf()[1]) == {"CustomGly", "CustomAla"}
        assert set(MonomerLibrary().GetMonomerNames()) == {"CustomGly", "CustomAla"}
        assert Sequence("CustomGly-CustomAla").is_valid()
    finally:
        monomer_store._invalidate_sdf()


def test_missing_library_override_fails_without_changing_explicit_registration(
    tmp_path,
    monkeypatch,
):
    from pyPept.interfaces.cli_monomer import register_monomer

    missing = tmp_path / "missing.sdf"
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(missing))
    with pytest.raises(ValueError, match="CABILN_MONOMER_LIBRARY.*existing"):
        Sequence("A")
    target = tmp_path / "explicit.sdf"
    register_monomer("NCC(=O)O", symbol="CustomGly", sdf_path=target)
    assert target.exists() and not missing.exists()


@pytest.mark.parametrize(
    "symbol", ["G", "A", "P", "H", "W", "K", "C", "ac", "am", "Bn_", "Mpa", "TATA"]
)
def test_standalone_monomer_display_matches_single_monomer_assembly(symbol):
    from pyPept.monomer_store import _load_sdf
    from pyPept.web.monomer_display import _restore_leaving_groups

    _, by_abbr = _load_sdf()
    restored = _restore_leaving_groups(by_abbr[symbol])
    assert Chem.MolToSmiles(restored) == smiles(symbol)
    assert all(atom.GetAtomicNum() != 0 for atom in restored.GetAtoms())


@pytest.mark.parametrize("leaving", ["[OH]", "[O-]", "[18OH]", "[18O-]"])
def test_display_and_assembly_preserve_terminal_charge_and_isotope(leaving):
    from pyPept.web.monomer_display import _restore_leaving_groups

    sequence = Sequence("A")
    sequence.s_monomers[0]["m_Rgroups"][1] = leaving
    template = Chem.Mol(sequence.s_monomers[0]["m_romol"])
    template.SetProp("m_Rgroups", f"[H],{leaving},[H]")
    expected = Chem.MolToSmiles(Chem.MolFromSmiles(f"N[C@@H](C)C(=O){leaving}"))
    assert Chem.MolToSmiles(Molecule(sequence).mol) == expected
    assert Chem.MolToSmiles(_restore_leaving_groups(template)) == expected


def test_display_and_assembly_keep_isotopic_hydrogen():
    from pyPept.web.monomer_display import _restore_leaving_groups

    sequence = Sequence("A")
    sequence.s_monomers[0]["m_Rgroups"][0] = "[2H]"
    template = Chem.Mol(sequence.s_monomers[0]["m_romol"])
    template.SetProp("m_Rgroups", "[2H],[OH],[H]")
    expected = Chem.MolToSmiles(Chem.MolFromSmiles("[2H]N[C@@H](C)C(=O)O"))
    assert Chem.MolToSmiles(Molecule(sequence).mol) == expected
    assert Chem.MolToSmiles(_restore_leaving_groups(template)) == expected


@pytest.mark.parametrize("leaving", ["OC", "notASmiles", "[F-]", "*"])
def test_display_rejects_invalid_leaving_groups_instead_of_truncating(leaving):
    from pyPept.web.monomer_display import _restore_leaving_groups

    sequence = Sequence("A")
    sequence.s_monomers[0]["m_Rgroups"][1] = leaving
    template = Chem.Mol(sequence.s_monomers[0]["m_romol"])
    template.SetProp("m_Rgroups", f"[H],{leaving},[H]")
    with pytest.raises(ValueError, match="leaving group|restoration|sanitization"):
        _restore_leaving_groups(template)
    with pytest.raises(ValueError, match="leaving group|restoration|sanitization"):
        Molecule(sequence)


def test_restoration_preserves_residue_tracking_and_pdb_atom_metadata():
    from pyPept.web.monomer_display import _restore_leaving_groups

    template = Chem.Mol(Sequence("A").s_monomers[0]["m_romol"])
    template.SetProp("m_Rgroups", "[H],[18O-],[H]")
    terminal = next(
        atom
        for atom in template.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 2
    )
    terminal.SetIntProp("_residue_idx", 7)
    terminal.SetMonomerInfo(
        Chem.AtomPDBResidueInfo(
            atomName=" OXT",
            residueName="ALA",
            residueNumber=8,
        )
    )
    restored = _restore_leaving_groups(template)
    oxygen = next(atom for atom in restored.GetAtoms() if atom.GetIsotope() == 18)
    assert oxygen.GetIntProp("_residue_idx") == 7
    assert oxygen.GetPDBResidueInfo().GetName() == " OXT"
    assert oxygen.GetPDBResidueInfo().GetResidueNumber() == 8
    assert terminal.GetAtomicNum() == 0 and terminal.GetIsotope() == 2


def test_reagent_preview_uses_shared_restoration_for_nonbonding_slots(monkeypatch):
    from pyPept.web import monomer_display

    template = Chem.Mol(Sequence("A").s_monomers[0]["m_romol"])
    template.SetProp("m_Rgroups", "[2H],[OH],[H]")
    monkeypatch.setattr(
        monomer_display,
        "_CAP_REACTIONS",
        {
            "CustomReagent": {"reagent_lg": "Cl", "reaction": "acylation"},
        },
    )
    reagent, metadata = monomer_display._restore_reagent_form(template, "CustomReagent")
    expected = Chem.MolToSmiles(Chem.MolFromSmiles("[2H]N[C@@H](C)C(=O)Cl"))
    assert Chem.MolToSmiles(reagent) == expected
    assert metadata["reagent_lg"] == "Cl"


@pytest.mark.parametrize("group", ["&1", "o1"])
@pytest.mark.parametrize("all_backbones", [False, True])
def test_activation_rejects_enhanced_stereo_groups(group, all_backbones):
    from pyPept.interfaces.monomer_pipeline import pre_activate_all

    activate = pre_activate_all if all_backbones else pre_activate
    with pytest.raises(ValueError, match="AND/OR stereochemistry"):
        activate(f"N[C@@H](C)C(=O)O |{group}:1|")


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_registration_rejects_enhanced_stereo_without_writing(tmp_path, group):
    from pyPept.monomer_store import register_molecule

    mol = Chem.MolFromSmiles(f"[1*]N([3*])[C@@H](C)C([2*])=O |{group}:3|")
    mol.SetProp("symbol", "RelativeAla")
    mol.SetProp("m_abbr", "RelativeAla")
    path = tmp_path / "custom.sdf"
    with pytest.raises(ValueError, match="AND/OR stereochemistry"):
        register_molecule(mol, sdf_path=path)
    assert not path.exists()


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_assembly_rejects_enhanced_stereo_templates(group):
    sequence = Sequence(f"<[1*]N[C@@H](C)C([2*])=O |{group}:2|>")
    with pytest.raises(ValueError, match="AND/OR stereochemistry"):
        Molecule(sequence)


def test_absolute_stereo_groups_remain_supported():
    result = pre_activate("N[C@@H](C)C(=O)O |a:1|")
    assert result.chuckles == pre_activate("N[C@@H](C)C(=O)O").chuckles
    assert smiles("<[1*]N[C@@H](C)C([2*])=O |a:2|>") == smiles("A")


@pytest.mark.parametrize("abbreviation", ["AliasGly", None])
def test_external_library_identifiers_resolve_consistently(
    tmp_path, monkeypatch, abbreviation
):
    from pyPept import monomer_store

    _, library = monomer_store._load_sdf()
    mol = Chem.Mol(library["G"])
    mol.SetProp("symbol", "CanonicalGly")
    if abbreviation is None:
        mol.ClearProp("m_abbr")
    else:
        mol.SetProp("m_abbr", abbreviation)
    path = tmp_path / "custom.sdf"
    path.write_text(Chem.SDWriter.GetText(mol, kekulize=False))
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    try:
        lookup = monomer_store._load_sdf()[1]
        expected_names = {"CanonicalGly", abbreviation or "CanonicalGly"}
        assert set(lookup) == expected_names
        expected = Chem.MolToSmiles(Chem.MolFromSmiles("NCC(=O)O"))
        for name in expected_names:
            assert smiles(name) == expected
        assert MonomerLibrary().GetMonomerNames() == (abbreviation or "CanonicalGly",)
    finally:
        monomer_store._invalidate_sdf()


@pytest.mark.parametrize(
    "identifiers",
    [
        [("One", "Shared"), ("Two", "Shared")],
        [("One", "Two"), ("Two", "Other")],
        [("Same", "AliasOne"), ("Same", "AliasTwo")],
    ],
)
def test_ambiguous_sdf_names_are_rejected(tmp_path, monkeypatch, identifiers):
    from pyPept import monomer_store

    _, library = monomer_store._load_sdf()
    records = []
    for source, (symbol, abbreviation) in zip(("G", "A"), identifiers):
        mol = Chem.Mol(library[source])
        mol.SetProp("symbol", symbol)
        mol.SetProp("m_abbr", abbreviation)
        records.append(Chem.SDWriter.GetText(mol, kekulize=False))
    path = tmp_path / "custom.sdf"
    path.write_text("".join(records))
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    try:
        with pytest.raises(ValueError, match="Ambiguous monomer name"):
            monomer_store._load_sdf()
        with pytest.raises(ValueError, match="Ambiguous monomer name"):
            Sequence(identifiers[0][0])
        with pytest.raises(ValueError, match="Ambiguous monomer name"):
            MonomerLibrary()
    finally:
        monomer_store._invalidate_sdf()


def test_ambiguous_csv_alias_cannot_silently_select_a_structure(tmp_path, monkeypatch):
    from pyPept.interfaces.cli_monomer import register_monomer

    path = tmp_path / "custom.sdf"
    register_monomer("NCC(=O)O", symbol="One", sdf_path=path)
    register_monomer("N[C@@H](C)C(=O)O", symbol="Two", sdf_path=path)
    (tmp_path / "monomers.csv").write_text(
        'token,synonyms\nOne,"Shared,OnlyOne,Two"\nTwo,Shared\n'
    )
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    assert smiles("OnlyOne") == Chem.MolToSmiles(Chem.MolFromSmiles("NCC(=O)O"))
    assert smiles("Two") == Chem.MolToSmiles(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    with pytest.raises(ValueError, match="Ambiguous monomer alias.*Shared"):
        Sequence("Shared")


@pytest.mark.parametrize("group", ["&1", "o1"])
def test_bulk_helm_import_skips_relative_stereo_and_keeps_absolute(tmp_path, group):
    import csv

    from pyPept.interfaces.monomer_pipeline import import_helm_sdf

    source, target = tmp_path / "source.sdf", tmp_path / "monomers.csv"
    records = []
    for symbol, stereo in (("RelativeAla", group), ("AbsoluteAla", "a")):
        mol = Chem.MolFromSmiles(f"N[C@@H](C)C(=O)O |{stereo}:1|")
        for name, value in {
            "polymerType": "PEPTIDE",
            "monomerType": "Backbone",
            "symbol": symbol,
            "name": symbol,
        }.items():
            mol.SetProp(name, value)
        records.append(Chem.SDWriter.GetText(mol, kekulize=False))
    source.write_text("".join(records))
    stored = list(Chem.SDMolSupplier(str(source)))
    assert stored[0].GetStereoGroups()[0].GetGroupType() in (
        Chem.StereoGroupType.STEREO_AND,
        Chem.StereoGroupType.STEREO_OR,
    )

    assert import_helm_sdf(source, target) == (1, 1)
    with target.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["token"] for row in rows] == ["AbsoluteAla"]
    restored = Chem.MolFromSmiles(rows[0]["input"])
    expected = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(expected)
