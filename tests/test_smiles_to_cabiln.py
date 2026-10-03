"""Reverse recognition preserves chemistry and editable atom partitions."""

import logging
import os
import sys

import pytest
from hypothesis import given, strategies as st
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.molecule import Molecule
from pyPept.sequence import Sequence

from _chemistry_oracles import _assert_same_monomer_partition
from _chemistry_fuzz import alternate_smiles, isolated_library, reordered, write_library
from _fuzzing import fuzz_settings, record


def _assert_fuzz_source_partition(source, result):
    """Check emitted ownership against independent source-to-product isomorphisms."""
    molecule = Chem.MolFromSmiles(source)
    assembly = Molecule(Sequence(result.cabiln), depiction=None)
    rebuilt = assembly.get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(molecule)
    owned = [atom for item in result.assignments for atom in item.source_atoms]
    assert sorted(owned) == list(range(molecule.GetNumAtoms()))
    groups = assembly.get_residue_atom_map()
    matches = rebuilt.GetSubstructMatches(molecule, useChirality=True, uniquify=False, maxMatches=4096)
    assert any(all(
        {match[atom] for atom in item.source_atoms} == set(groups[item.residue_index])
        for item in result.assignments
    ) for match in matches)
    for item in result.assignments:
        assert all(atom in item.source_atoms for slot, atom in item.attachments)


@pytest.mark.fuzz
@fuzz_settings(examples=10)
@given(glycines=st.integers(2, 3), length=st.integers(1, 3), aliases=st.integers(0, 2),
       overlap=st.booleans(), order=st.integers(0, 65535), notation=st.sampled_from(("percent", "bracket")))
def test_fuzz_known_unknown_partitions_survive_atom_library_and_alias_order(
    glycines, length, aliases, overlap, order, notation
):
    from pyPept.interfaces.cli_monomer import register_monomer
    from pyPept.recognition import RecognitionBudgets
    from pyPept.smiles import convert_smiles

    original = "NCC(=O)" * glycines + f"N[C@@H]({'C' * length}(F)(F)F)C(=O)N[C@@H](C)C(=O)O"
    source = alternate_smiles(original, order)
    record("recognition.known-unknown", source=source, aliases=aliases, overlap=overlap,
           order=order, notation=notation)
    with isolated_library(("G", "A")) as path:
        records = list(Chem.SDMolSupplier(str(path), removeHs=False))
        glycine = next(item for item in records if item.GetProp("m_abbr") == "G")
        for index in range(aliases):
            alias = Chem.Mol(glycine)
            alias.SetProp("symbol", f"GlyAlias{index}")
            alias.SetProp("m_abbr", f"GlyAlias{index}")
            records.append(alias)
        write_library(path, records)
        if overlap:
            register_monomer("NCC(=O)NCC(=O)O", "GlyPair")
        records = list(Chem.SDMolSupplier(str(path), removeHs=False))
        write_library(path, (reordered(item, order) for item in reversed(records)))
        result = convert_smiles(source, notation=notation, budgets=RecognitionBudgets(max_states=5000))
        assert result.recognition_status == "partial"
        assert not result.inferred_stereo
        _assert_fuzz_source_partition(source, result)
        molecule = Chem.MolFromSmiles(source)
        original_molecule = Chem.MolFromSmiles(original)
        match = molecule.GetSubstructMatch(original_molecule, useChirality=True)
        assert len(match) == molecule.GetNumAtoms()
        expected_known = {match[index] for index in range(4 * glycines)} | set(match[-6:])
        known = {atom for item in result.assignments if item.recognized for atom in item.source_atoms}
        assert known == expected_known
        assert sum(not item.recognized for item in result.assignments) == 1
        assert any(item.recognized and item.symbol == "A" for item in result.assignments)
        assert result.synthetic_components == (0,)


@pytest.mark.fuzz
@pytest.mark.parametrize("nitrogen,carbon", [("[NH3+]", "C"), ("N", "[13CH2]"), ("[NH3+]", "[13CH2]")])
@fuzz_settings(examples=10)
@given(order=st.integers(0, 65535), notation=st.sampled_from(("percent", "bracket")))
def test_fuzz_charge_and_isotope_changes_are_preserved_and_not_claimed_as_glycine(
    nitrogen, carbon, order, notation
):
    from pyPept.smiles import convert_smiles

    source = alternate_smiles(f"{nitrogen}{carbon}C(=O)N[C@@H](C)C(=O)O", order)
    record("recognition.charge-isotope", source=source, notation=notation)
    with isolated_library(("G", "A")):
        result = convert_smiles(source, notation=notation)
        assert result.recognition_status == "partial"
        _assert_fuzz_source_partition(source, result)
        assert [item.symbol for item in result.assignments if item.recognized] == ["A"]
        assert sum(not item.recognized for item in result.assignments) == 1


@pytest.mark.fuzz
@pytest.mark.parametrize("stereo", ("@", "@@", ""))
@fuzz_settings(examples=10)
@given(order=st.integers(0, 65535), allow_inference=st.booleans())
def test_fuzz_stereo_inference_is_explicit_and_never_overrides_defined_source(
    stereo, order, allow_inference
):
    from pyPept.recognition import RecognitionBudgets
    from pyPept.smiles import convert_smiles

    source = alternate_smiles(f"N[C{stereo}H](C)C(=O)O" if stereo else "NC(C)C(=O)O", order)
    record("recognition.stereo-policy", source=source, allow_inference=allow_inference)
    with isolated_library(("A",)):
        result = convert_smiles(source, budgets=RecognitionBudgets(allow_unspecified_stereo=allow_inference))
        inferred = not stereo and allow_inference
        assert result.inferred_stereo is inferred
        assert any("stereochemistry unspecified" in warning for warning in result.warnings) is inferred
        known = stereo == "@@" or inferred
        assert any(item.recognized for item in result.assignments) is known
        if inferred:
            assert Chem.MolToSmiles(Molecule(Sequence(result.cabiln), depiction=None).get_molecule(fmt="ROMol")) == Chem.MolToSmiles(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
        else:
            _assert_fuzz_source_partition(source, result)


@pytest.mark.fuzz
@fuzz_settings(examples=10)
@given(length=st.integers(2, 5), order=st.integers(0, 65535))
def test_fuzz_exhausted_recognition_preserves_structure_without_claiming_completion(length, order):
    from pyPept.recognition import RecognitionBudgets
    from pyPept.smiles import convert_smiles

    source = alternate_smiles("NCC(=O)" * length + "O", order)
    record("recognition.budget", source=source, max_states=1)
    with isolated_library(("G",)):
        result = convert_smiles(source, budgets=RecognitionBudgets(max_states=1))
        assert result.search_complete is False
        assert any("limit" in warning.lower() or "budget" in warning.lower() for warning in result.warnings)
        rebuilt = Molecule(Sequence(result.cabiln), depiction=None).get_molecule(fmt="ROMol")
        assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(Chem.MolFromSmiles(source))


@pytest.mark.fuzz
@pytest.mark.parametrize("symbols,product", [
    (("Ala_3Br", "C"), "N[C@@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    (("D_Ala_3Cl", "C"), "N[C@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    (("C", "C"), "N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O"),
    (("Pra", "AzK"), "N[C@@H](Cc1cn(CCCC[C@H](N)C(=O)O)nn1)C(=O)O"),
])
@fuzz_settings(examples=10)
@given(order=st.integers(0, 65535), alias=st.booleans(),
       notation=st.sampled_from(("percent", "bracket")))
def test_fuzz_non_amide_recognition_retains_editable_partitions(symbols, product, order, alias, notation):
    from pyPept.smiles import convert_smiles

    source = alternate_smiles(product, order)
    record("recognition.non-amide", symbols=symbols, source=source, alias=alias, notation=notation)
    with isolated_library(tuple(dict.fromkeys(symbols))) as path:
        records = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if alias:
            extra = Chem.Mol(records[0])
            extra.SetProp("symbol", "ZZAlias")
            extra.SetProp("m_abbr", "ZZAlias")
            records.append(extra)
        write_library(path, (reordered(item, order) for item in reversed(records)))
        result = convert_smiles(source, notation=notation)
        assert result.recognition_status == "complete"
        assert len(result.assignments) == 2
        assert all(item.recognized for item in result.assignments)
        assert not result.synthetic_components
        _assert_fuzz_source_partition(source, result)
        _assert_same_monomer_partition(
            f"{symbols[0]}.!bond(4,4)%{symbols[1]}.!bond(4,4)", result.cabiln)


@pytest.mark.parametrize("template,leaving,diagnostic", [
    pytest.param("[1*]C[1*]", "[H]", "attachment slot numbers must be unique",
                 id="invalid-template"),
    pytest.param("[1*]C", "[Cl]C", "slot mask 0", id="invalid-template-state"),
])
def test_skipped_library_template_diagnostics_are_not_import_warnings(
    monkeypatch, caplog, template, leaving, diagnostic,
):
    from pyPept import recognition
    from pyPept.smiles import convert_smiles

    stamp, molecules = recognition._snapshot()
    broken = Chem.MolFromSmiles(template)
    broken.SetProp("symbol", "InvalidUnusedTemplate")
    broken.SetProp("m_Rgroups", leaving)
    monkeypatch.setattr(
        recognition, "_snapshot",
        lambda: ((stamp, "invalid-unused-template", template), (*molecules, broken)),
    )

    with caplog.at_level(logging.DEBUG, logger="pyPept.recognition"):
        result = convert_smiles("NCC(=O)O")

    assert result.cabiln == "G"
    assert result.recognition_status == "complete"
    assert result.search_complete
    assert result.warnings == ()
    assert "InvalidUnusedTemplate" in caplog.text
    assert diagnostic in caplog.text


@pytest.mark.parametrize("symbol,source", [
    ("Ala_3Br", "N[C@@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    ("Ala_3Cl", "N[C@@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    ("ClAcAla", "N[C@@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    ("D_Ala_3Br", "N[C@H](CSC[C@H](N)C(=O)O)C(=O)O"),
    ("D_Ala_3Cl", "N[C@H](CSC[C@H](N)C(=O)O)C(=O)O"),
])
def test_corrected_halide_substitution_is_recognized_with_exact_ownership(
    tmp_path, monkeypatch, symbol, source,
):
    from pyPept.monomer_store import _load_sdf
    from pyPept.smiles import convert_smiles

    library = _load_sdf()[1]
    path = tmp_path / "halide-cysteine.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for name in (symbol, "C"):
            writer.write(library[name])
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    result = convert_smiles(source)
    assert result.recognition_status == "complete"
    assert result.synthetic_components == ()
    assert {item.symbol for item in result.assignments} == {symbol, "C"}
    rebuilt = Molecule(Sequence(result.cabiln)).get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(Chem.MolFromSmiles(source))
    _assert_same_monomer_partition(f"{symbol}.!s(4,4)%C.!s(4,4)", result.cabiln)


def test_partial_peptide_does_not_exhaust_search_on_unsupported_boundaries(monkeypatch):
    from pyPept import recognition
    from pyPept.smiles import convert_smiles

    # A source-specific histidine tautomer plus eight ordinary glutamates. The
    # complete library includes haloalanine cores that also match pieces of E;
    # their unsupported C-C connections must not swamp the useful partition.
    source = (
        "N[C@@H](Cc1c[nH]cn1)C(=O)"
        + "N[C@@H](CCC(=O)O)C(=O)" * 8 + "O"
    )
    original_search = recognition._search
    states = []

    def search(*args, **kwargs):
        result = original_search(*args, **kwargs)
        states.append(result[2])
        return result

    monkeypatch.setattr(recognition, "_search", search)
    result = convert_smiles(source, budgets=recognition.RecognitionBudgets(max_states=5000))
    assert [symbol for symbol, _, _ in result.details[1:]] == ["E"] * 8
    assert result.details[0][0].startswith("<")
    assert result.synthetic_components == (0,)
    rebuilt = Molecule(Sequence(result.cabiln)).get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(Chem.MolFromSmiles(source))
    assert states and max(states) < 5000


def _s2c_roundtrip(biln: str):
    """CABILN -> SMILES -> CABILN.

    Returns (cabiln_str, details) where details is the list of
    (abbr, n_matched, n_total_atoms) tuples from smiles_to_cabiln_core.
    """
    from pyPept.smiles import smiles_to_cabiln_core  # noqa: E402
    mol = Molecule(Sequence(biln)).get_molecule(fmt='ROMol')
    smi = Chem.MolToSmiles(mol)
    result, details = smiles_to_cabiln_core(smi)
    _assert_same_monomer_partition(biln, result)
    return result, details


class TestSmilesToCabiln:
    """SMILES conversion preserves chemistry and editable monomer boundaries.

    Literal notation is checked for simple, unambiguous cases. Branches and
    symmetric scaffolds also require the same atom partition after assembly.
    """

    @staticmethod
    def _r(biln: str):
        """Shorthand that returns (cabiln, details)."""
        return _s2c_roundtrip(biln)

    # ------------------------------------------------------------------
    # Terminal caps and residue boundaries
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("biln,expected_abbr,expected_cabiln", [
        # Single residue, C-cap (am) only.
        pytest.param("G-am", "G", "G-am", id="single_G_am_only"),
        pytest.param("A-am", "A", "A-am", id="single_A_am_only"),
        pytest.param("V-am", "V", "V-am", id="single_V_am_only"),

        # Single residue, N-cap (ac) only.
        pytest.param("ac-G", "G", "ac-G", id="single_ac_G_only"),
        pytest.param("ac-A", "A", "ac-A", id="single_ac_A_only"),

        # Both caps must remain separate from the single backbone occurrence.
        pytest.param("ac-G-am", "G", "ac-G-am", id="single_ac_G_am_both"),
        pytest.param("ac-A-am", "A", "ac-A-am", id="single_ac_A_am_both"),
        pytest.param("ac-V-am", "V", "ac-V-am", id="single_ac_V_am_both"),

        # Multi-residue chains retain both caps and each backbone occurrence.
        pytest.param("ac-A-G-am",     "A",  "ac-A-G-am",     id="multi_both_caps_2mer"),
        pytest.param("ac-A-G-V-L-am", "A",  "ac-A-G-V-L-am", id="multi_both_caps_4mer"),

        # Multi-residue, am cap only.
        pytest.param("A-G-am",   "A",  "A-G-am",   id="multi_am_only_2mer"),
        pytest.param("A-G-V-am", "A",  "A-G-V-am", id="multi_am_only_3mer"),

        # Multi-residue, ac cap only.
        pytest.param("ac-A-G",   "A",  "ac-A-G",   id="multi_ac_only_2mer"),
        pytest.param("ac-A-G-V", "A",  "ac-A-G-V", id="multi_ac_only_3mer"),

        # No caps: output equals input, no '?' residues.
        pytest.param("A-G", "A", "A-G", id="no_caps_2mer"),
        pytest.param("G",   "G", "G",   id="no_caps_single"),
    ])
    def test_cap_stripping(self, biln, expected_abbr, expected_cabiln):
        """N- and C-caps stay separate from each residue, including a single residue."""
        result, details = self._r(biln)
        assert result == expected_cabiln, (
            f"Full CABILN mismatch for {biln!r}:\n"
            f"  expected: {expected_cabiln!r}\n"
            f"  got:      {result!r}"
        )
        first_abbr = details[0][0]
        assert first_abbr == expected_abbr, (
            f"First-residue abbr wrong for {biln!r}: "
            f"expected {expected_abbr!r}, got {first_abbr!r} (details={details})"
        )
        assert "?" not in result, f"Unknown residue '?' in: {result!r}"


    # ------------------------------------------------------------------
    # Library residue identity and specified stereochemistry
    # ------------------------------------------------------------------

    def test_all_20_standard_aa(self):
        """All 20 standard amino acids retain their canonical library symbols."""
        biln = "A-G-V-L-I-P-F-W-M-S-T-C-Y-H-K-R-D-E-N-Q-am"
        result, details = self._r(biln)
        # The output format is "A-G-V-…-Q-am"; split and drop the trailing am.
        tokens = result.split("-")
        # Last token is 'am'; first 20 are the backbone residues.
        assert tokens[-1] == "am", f"Expected trailing 'am', got: {tokens[-1]!r}"
        backbone = tokens[:20]
        expected = ["A","G","V","L","I","P","F","W","M","S",
                    "T","C","Y","H","K","R","D","E","N","Q"]
        assert backbone == expected, (
            f"Canonical AA mismatch:\n  expected: {expected}\n  got: {backbone}"
        )


    def test_d_amino_acid_stereo_disambiguation(self):
        """Specified D stereochemistry selects D-amino-acid library templates."""
        result, details = self._r("dA-dV-G-am")
        assert "?" not in result, f"Unknown residue in: {result!r}"
        abbrs = [d[0] for d in details]
        # Library stores D-Ala as 'DAla', D-Val as 'DVal'.
        assert abbrs[0] in {"dA", "DAla"}, (
            f"First residue should be D-Ala variant, got {abbrs[0]!r}"
        )
        assert abbrs[1] in {"dV", "DVal"}, (
            f"Second residue should be D-Val variant, got {abbrs[1]!r}"
        )
        assert abbrs[2] == "G", f"Third residue should be G, got {abbrs[2]!r}"

    def test_non_natural_backbone_monomers(self):
        """Aib and Orn remain recognizable backbone monomer occurrences."""
        result, details = self._r("Aib-Orn-am")
        assert "?" not in result, f"Unknown residue in: {result!r}"
        abbrs = [d[0] for d in details]
        assert abbrs[0] == "Aib", f"Expected Aib, got {abbrs[0]!r}"
        assert abbrs[1] == "Orn", f"Expected Orn, got {abbrs[1]!r}"

    def test_l_vs_d_leucine_stereo(self):
        """L- and D-leucine have distinct recognized symbols and exact stereo."""
        result_l,  details_l  = self._r("L-am")
        result_dl, details_dl = self._r("dL-am")
        assert "?" not in result_l,  f"Unknown in L-am: {result_l!r}"
        assert "?" not in result_dl, f"Unknown in dL-am: {result_dl!r}"
        assert details_l[0][0] != details_dl[0][0], (
            f"L-Leu and D-Leu produced same abbr: {details_l[0][0]!r}"
        )

    # ------------------------------------------------------------------
    # Cyclic peptide detection
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("biln", [
        # Minimal 3-mer cyclic (head-to-tail backbone ring).
        pytest.param("!1-A-G-K-!1",   id="cyclic_3mer"),
        # 4-mer cyclic.
        pytest.param("!1-A-G-K-D-!1", id="cyclic_4mer"),
    ])
    def test_cyclic_detection(self, biln):
        """A head-to-tail cycle is emitted with paired terminal ring markers."""
        result, details = self._r(biln)
        assert result.startswith("!1-"), f"Missing cyclic prefix: {result!r}"
        assert result.endswith("-!1"),   f"Missing cyclic suffix: {result!r}"
        assert "?" not in result,        f"Unknown residue in cyclic: {result!r}"

    # ------------------------------------------------------------------
    # Branched peptides and lipid linkers
    # ------------------------------------------------------------------

    def test_minimal_lipid_linker_bracket(self):
        """A lipid-bearing Lys preserves the full branch and its flanking residues."""
        biln   = "A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in lipid-linker result: {result!r}"
        assert "K.[" in result or "K.!1" in result, (
            f"Expected bracket/crosslink notation on K in: {result!r}"
        )

    def test_retatrutide_full_sequence(self):
        """The 39-residue model retains its backbone and complete C20 lipid branch."""
        biln = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
            "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
        )
        result, details = self._r(biln)

    # ------------------------------------------------------------------
    # Crosslink and staple annotation
    # ------------------------------------------------------------------

    def test_staple_i_i4_same_handedness(self):
        """An i,i+4 S5/S5 metathesis staple retains both R4 endpoints and terminal caps."""
        biln   = "ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in staple result: {result!r}"
        assert result.count(".!1") == 2, (
            f"Expected exactly two .!1 annotations in: {result!r}"
        )
        assert "(4,4)" in result, f"Missing (4,4) r-group annotation in: {result!r}"
        assert result.startswith("ac-"), f"N-cap lost: {result!r}"
        assert result.endswith("-am"),   f"C-cap lost: {result!r}"

    def test_staple_i_i7_opposite_handedness(self):
        """An i,i+7 S5/R8 staple retains both endpoints and their opposite configurations."""
        biln   = "ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in i+7 staple result: {result!r}"
        assert result.count(".!1") == 2, (
            f"Expected two .!1 annotations in: {result!r}"
        )

    # ------------------------------------------------------------------
    # Combined topology features
    # ------------------------------------------------------------------

    def test_priority_ac_am_lipid_branch(self):
        """Both terminal caps, the A/K/G backbone, and the complete lipid branch survive."""
        biln   = "ac-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue: {result!r}"
        assert "K.[" in result or "K.!1" in result, (
            f"Lipid branch on K lost: {result!r}"
        )


    def test_priority_cyclic_plus_lipid_branch(self):
        """A cyclic backbone and its complete lipid branch survive together."""
        biln   = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in cyclic+branch: {result!r}"
        assert result.startswith("!1-"), f"Cyclic prefix lost: {result!r}"
        assert result.split("%")[0].endswith("-!1"),   f"Cyclic suffix lost: {result!r}"
        assert "K.[" in result or "K.!" in result, (
            f"Lipid branch on K lost in cyclic+branch: {result!r}"
        )

    def test_dual_lipid_branch_terminal_k(self):
        """Both terminal Lys occurrences retain their separate E/AEEA/C20FA arms."""
        biln = "ac-K.[E(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E(4,4).AEEA(1,2).C20FA(1,2)]-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in dual-lipid result: {result!r}"
        assert result.count("K.!") == 2, (
            f"Expected two attached K residues, got: {result!r}"
        )

    # ------------------------------------------------------------------
    # TBMB three-way thioether crosslinker
    # ------------------------------------------------------------------

    def test_tbmb_assembly_three_cys(self):
        """TBMB scaffold: assembly of three Cys thioether crosslinks succeeds.

        ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3
        All three thioether bonds (Cys-S → TBMB-CH2) must form, yielding a
        closed tricyclic structure with 51 heavy atoms.
        """
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        assert romol is not None, "TBMB assembly returned None"
        assert romol.GetNumAtoms() == 51, (
            f"Expected 51 heavy atoms, got {romol.GetNumAtoms()}"
        )
        from rdkit import Chem
        smi = Chem.MolToSmiles(romol)
        assert smi, "TBMB assembly produced empty SMILES"
        # Verify three thioether S-C bonds are present (3× Cys-S-CH2-Ar)
        assert smi.count('CSC') == 3, (
            f"Expected 3 thioether CSC motifs in TBMB product, got {smi.count('CSC')}: {smi}"
        )

    def test_tbmb_converter_roundtrip(self):
        """TBMB: branch -> bracket (flat) -> branch roundtrip, and flat bracket -> branch."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        flat   = "ac-C.[TBMB(4,4).!2(5,4).!3(6,4)]-A-A-C.!2-A-A-C.!3-am"
        bracket = cabiln_to_bracket(branch)
        assert bracket == flat, f"branch->bracket: {bracket!r}"
        back = cabiln_to_branch(bracket)
        assert back == branch, f"bracket->branch: {back!r}"

    def test_tbmb_assembly_partial_two_cys(self):
        """TBMB scaffold: partial use (two of three Cys crosslinks) also assembles.

        ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2 leaves one CH2Br unreacted.
        """
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2"
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        assert romol is not None, "Partial TBMB assembly returned None"
        from rdkit import Chem
        smi = Chem.MolToSmiles(romol)
        assert smi, "Partial TBMB assembly produced empty SMILES"
        assert smi.count('CSC') == 2, (
            f"Expected 2 thioether CSC motifs, got {smi.count('CSC')}: {smi}"
        )
        assert 'CBr' in smi, f"Expected unreacted CH2Br (CBr) in partial TBMB: {smi}"

    def test_tbmb_smiles_to_cabiln_roundtrip(self):
        """TBMB SMILES→CABILN round-trip."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        result, _ = self._r(cabiln)

    def test_tbmb_asymmetric_roundtrip(self):
        """Asymmetric-loop TBMB bicycle SMILES→CABILN round-trip."""
        cabiln = "ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TBMB.!1.!2.!3"
        result, _ = self._r(cabiln)


    def test_tbmb_scaffold_plus_disulfide_roundtrip(self):
        """TBMB scaffold + independent Cys-Cys disulfide on same peptide."""
        cabiln = "ac-C.!1(4,4)-A-C.!4(4,4)-A-C.!4(4,4)-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        result, _ = self._r(cabiln)

    # ── TATA scaffold (thio-Michael) ──────────────────────────────────────────

    def test_tata_assembly_three_cys(self):
        """TATA scaffold: assembly of three Cys thio-Michael crosslinks succeeds."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TATA.!1.!2.!3"
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        assert romol is not None
        smi = Chem.MolToSmiles(romol)
        assert '.' not in smi, "TATA product is disconnected"
        assert smi.count('CSC') == 3, f"Expected 3 CSC thioether motifs, got {smi.count('CSC')}"

    def test_tata_smiles_to_cabiln_roundtrip(self):
        """TATA SMILES→CABILN round-trip."""
        cabiln = "ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TATA.!1.!2.!3"
        result, _ = self._r(cabiln)

    # ── TBAB scaffold (benzene triamide, alkyl halide arms) ───────────────────

    def test_tbab_assembly_three_cys(self):
        """TBAB scaffold: assembly of three Cys thioether crosslinks succeeds."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBAB.!1.!2.!3"
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        assert romol is not None
        smi = Chem.MolToSmiles(romol)
        assert '.' not in smi, "TBAB product is disconnected"
        assert smi.count('CSC') == 3, f"Expected 3 CSC thioether motifs, got {smi.count('CSC')}"

    def test_tbab_smiles_to_cabiln_roundtrip(self):
        """TBAB SMILES→CABILN round-trip."""
        cabiln = "ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TBAB.!1.!2.!3"
        result, _ = self._r(cabiln)

    # ── Gap 1: TATA 2-arm partial roundtrip ───────────────────────────────────

    def test_tata_partial_two_arm_roundtrip(self):
        """A TATA scaffold with two reacted arms retains its product and monomer boundaries."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TATA.!1.!2"
        result, _ = self._r(cabiln)

    # ── Gap 2: TBAB 2-arm partial roundtrip ───────────────────────────────────

    def test_tbab_partial_two_arm_roundtrip(self):
        """A TBAB scaffold with two reacted arms retains its product and monomer boundaries."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBAB.!1.!2"
        result, _ = self._r(cabiln)

    # ── Gap 3: Bracket notation (→ []) assembles same SMILES ──────────────────

    @pytest.mark.parametrize("pct_cabiln", [
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3", id="tbmb_3arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2",                   id="tbmb_2arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!1(4,4)-am",                               id="disulfide"),
    ])
    def test_bracket_notation_assembles_same_smiles(self, pct_cabiln):
        """→ [] path: bracket CABILN assembles to the identical SMILES as the percent form."""
        from pyPept.inputs import format_source

        ref_mol = Molecule(Sequence(pct_cabiln)).get_molecule(fmt='ROMol')
        ref_smi = Chem.MolToSmiles(ref_mol)

        bracket_cabiln = format_source(pct_cabiln, 'bracket')
        brk_mol = Molecule(Sequence(bracket_cabiln)).get_molecule(fmt='ROMol')
        brk_smi = Chem.MolToSmiles(brk_mol)

        assert brk_smi == ref_smi, (
            f"Bracket form yields different SMILES for {pct_cabiln!r}:\n"
            f"  percent:  {ref_smi}\n"
            f"  bracket:  {brk_smi}\n"
            f"  (bracket CABILN: {bracket_cabiln!r})"
        )

    # ── Gap 3b: Bracket form round-trips back to percent via SMILES→CABILN ───

    @pytest.mark.parametrize("pct_cabiln", [
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3", id="tbmb_3arm_rtrip"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2",                   id="tbmb_2arm_rtrip"),
    ])
    def test_bracket_to_percent_roundtrip(self, pct_cabiln):
        """Bracket CABILN assembled to SMILES then decoded via SMILES→CABILN gives percent form."""
        from pyPept.inputs import format_source
        from pyPept.smiles import smiles_to_cabiln_core

        bracket_cabiln = format_source(pct_cabiln, 'bracket')
        brk_mol = Molecule(Sequence(bracket_cabiln)).get_molecule(fmt='ROMol')
        brk_smi = Chem.MolToSmiles(brk_mol)
        result, _ = smiles_to_cabiln_core(brk_smi)

        _assert_same_monomer_partition(pct_cabiln, result)

    # ── Gap 4: CuAAC assembly produces a triazole ring ────────────────────────

    def test_cuaac_assembly_forms_triazole(self):
        """CuAAC: Pra (terminal alkyne R4) + AzK (azide R4) forms a 1,2,3-triazole."""
        cabiln = "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am"
        seq = Sequence(cabiln)
        mol = Molecule(seq).get_molecule(fmt='ROMol')
        assert mol is not None, "CuAAC assembly returned None"
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi, f"CuAAC product is disconnected: {smi!r}"
        triazole_smarts = Chem.MolFromSmarts('[n]1[n][n][c][c]1')
        assert mol.HasSubstructMatch(triazole_smarts), (
            f"No 1,2,3-triazole ring in CuAAC product: {smi!r}"
        )

    def test_cuaac_smiles_to_cabiln_roundtrip(self):
        """CuAAC product recognition retains the precursors and their reaction closure."""
        cabiln = "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am"
        result, _ = self._r(cabiln)
        assert result == cabiln

    # ── Gap 4b: Backbone-less fallback — cap-only and standalone monomers ─────

    @pytest.mark.parametrize("cabiln", [
        pytest.param("ac-am",   id="ac_am"),
        pytest.param("fmoc-am", id="fmoc_am"),
        pytest.param("boc-am",  id="boc_am"),
    ])
    def test_caponly_chain_roundtrip(self, cabiln):
        """Two-cap chains (no backbone residues) round-trip through SMILES."""
        result, _ = self._r(cabiln)

    def test_standalone_tbmb_smiles_detected(self):
        """Unreacted TBMB scaffold (all three Br arms present) is identified as 'TBMB'."""
        from pyPept.smiles import smiles_to_cabiln_core
        # Unreacted TBMB: three CH2Br arms, no Cys thioether bonds formed
        smi = 'BrCc1cc(CBr)cc(CBr)c1'
        result, _ = smiles_to_cabiln_core(smi)
        assert result == 'TBMB', f"Expected 'TBMB', got {result!r}"

    # ── Gap 5: Cyclic backbone + sidechain disulfide ──────────────────────────

    def test_cyclic_backbone_with_disulfide_roundtrip(self):
        """Head-to-tail cyclic peptide with an internal Cys-Cys disulfide round-trips.

        The ring-opening point and crosslink ID numbering may differ from the
        input (canonical form), so we verify structural invariants rather than
        exact string equality.
        """
        import re as _re
        cabiln = "!1-A-C.!2(4,4)-G-C.!2(4,4)-A-!1"
        result, _ = self._r(cabiln)
        assert result.startswith("!1-"),  f"Cyclic prefix lost: {result!r}"
        assert result.endswith("-!1"),    f"Cyclic suffix lost: {result!r}"
        assert "?" not in result,         f"Unknown residue in cyclic+disulfide: {result!r}"
        # Both Cys residues must carry a disulfide annotation (4,4); ID may be renumbered
        disulfide_hits = _re.findall(r'\.!\d+\(4,4\)', result)
        assert len(disulfide_hits) == 2, (
            f"Expected 2 Cys-Cys disulfide annotations in cyclic+disulfide: {result!r}"
        )

    # ── Gap 6: Stereo (L vs D) preserved through SMILES→CABILN ───────────────

    def test_l_and_d_stereo_preserved(self):
        """Specified tetrahedral configurations survive the complete molecular round trip."""
        for cabiln in ("ac-A-V-L-I-am", "ac-dA-G-V-am"):
            mol_before = Molecule(Sequence(cabiln)).get_molecule(fmt='ROMol')
            smi_before = Chem.MolToSmiles(mol_before)
            assert '@' in smi_before, f"No stereocenters in {cabiln!r} SMILES: {smi_before!r}"
            result, _ = self._r(cabiln)
            mol_after = Molecule(Sequence(result)).get_molecule(fmt='ROMol')
            smi_after = Chem.MolToSmiles(mol_after)
            assert smi_before == smi_after, (
                f"Stereo changed in round-trip for {cabiln!r}:\n"
                f"  before: {smi_before}\n"
                f"  after:  {smi_after}"
            )

    # ── Gap 7: Full scaffold match preferred over partial ─────────────────────

    def test_tata_full_match_preferred_over_partial(self):
        """3-arm TATA product: full (3-arm) scaffold detected, not a 2-arm partial."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TATA.!1.!2.!3"
        result, _ = self._r(cabiln)
        assert result.split("%TATA.")[1].count("!") == 3, f"3-arm TATA suffix missing: {result!r}"

    # ── Gap 8: thia_michael_c SMARTS specificity ─────────────────────────────

    def test_thia_michael_c_smarts_specificity(self):
        """thia_michael_c pre_smarts matches beta-acrylamide C but not thiol, carboxyl, or alkene."""
        from pyPept.site_chemistry import SITE_RULES

        entry = next(rule for rule in SITE_RULES if rule.name == 'thia_michael_c')
        thia_smarts = Chem.MolFromSmarts(entry.raw_smarts)
        assert thia_smarts is not None, "thia_michael_c pre_smarts failed to compile"

        for smi_match in ('C=CC(=O)N', 'C=CC(=O)NC', 'C=CC(=O)N1CCCC1'):
            mol = Chem.MolFromSmiles(smi_match)
            assert mol is not None and mol.HasSubstructMatch(thia_smarts), (
                f"thia_michael_c SMARTS should match acrylamide {smi_match!r}"
            )
        for smi_no in ('NCC(S)C(=O)O', 'NCC(CC(=O)O)C(=O)O', 'C=C', 'C=CC=O', 'C=CC(=O)O'):
            mol = Chem.MolFromSmiles(smi_no)
            assert mol is not None and not mol.HasSubstructMatch(thia_smarts), (
                f"thia_michael_c SMARTS must NOT match {smi_no!r}"
            )

    # ── Gap 9 & 10: Partial scaffold assembly produces no dummy atoms ─────────

    @pytest.mark.parametrize("pct_cabiln", [
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TATA.!1.!2", id="tata_2arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2", id="tbmb_2arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBAB.!1.!2", id="tbab_2arm"),
    ])
    def test_partial_scaffold_no_dummy_atoms(self, pct_cabiln):
        """2-of-3 scaffold arms reacted: assembled SMILES has no unresolved dummy (*) and is connected."""
        seq = Sequence(pct_cabiln)
        mol = Molecule(seq).get_molecule(fmt='ROMol')
        assert mol is not None, f"Assembly returned None for {pct_cabiln!r}"
        smi = Chem.MolToSmiles(mol)
        assert '*' not in smi, f"Dummy atom (*) in 2-arm partial product: {smi!r}"
        assert '.' not in smi, f"Disconnected 2-arm partial product: {smi!r}"

    # ── Cyclic backbone + sidechain crosslink: !n ID collision fix ────────────

    def test_cyclic_backbone_plus_disulfide_no_id_collision(self):
        """Cyclic backbone ring uses !1; sidechain crosslinks must start from !2.

        Before the fix, smiles_to_cabiln_core assigned !1 to both the ring
        closure and the disulfide, producing an un-reassemblable CABILN.
        After the fix the sidechain crosslinks are numbered from !2 onwards.
        """
        import re as _re
        input_biln = '!1-A-C.!2(4,4)-G-C.!2(4,4)-A-!1'
        result, _ = self._r(input_biln)
        # The ring-closure markers must be the only bare !1 at boundaries
        assert result.startswith('!1-') and result.endswith('-!1'), (
            f"Ring closure markers lost: {result!r}"
        )
        # Sidechain crosslinks must NOT use !1 (reserved for ring closure)
        sidechain_ids = _re.findall(r'\.!(\d+)\(', result)
        assert '1' not in sidechain_ids, (
            f"Sidechain crosslink incorrectly reused !1 (ring ID): {result!r}"
        )
        # Must reassemble cleanly
        mol2 = Molecule(Sequence(result)).get_molecule(fmt='ROMol')
        assert mol2 is not None, f"Reassembly of {result!r} returned None"
        smi2 = Chem.MolToSmiles(mol2)
        assert '.' not in smi2, f"Reassembled molecule is disconnected: {smi2!r}"

    # ------------------------------------------------------------------
    # N-cap identification (all types, not just 'ac')
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("biln,expected_cap_token", [
        pytest.param("ac-A-G-am",   "ac",   id="n_cap_ac"),
        pytest.param("boc-A-G-am",  "boc",  id="n_cap_boc"),
        pytest.param("fmoc-A-G-am", "fmoc", id="n_cap_fmoc"),
        pytest.param("Cbz-A-G-am",  "Cbz",  id="n_cap_Cbz"),
    ])
    def test_n_cap_type_identification(self, biln, expected_cap_token):
        """Each terminal protecting group is recognized from the current library."""
        result, details = self._r(biln)
        assert result.split("-")[0].casefold() == expected_cap_token.casefold(), (
            f"Expected N-cap {expected_cap_token!r} in output, got: {result!r}"
        )
        assert "?" not in result, f"Unknown residue in capped peptide: {result!r}"

    @pytest.mark.parametrize(
        'biln',
        [
            pytest.param('ac-A-NHMe', id='c_cap_NHMe'),
            pytest.param('ac-A-NHEt', id='c_cap_NHEt'),
            pytest.param('ac-A-OEt', id='c_cap_OEt'),
            pytest.param('ac-A-OtBu', id='c_cap_OtBu'),
            pytest.param('ac-A-_OMe', id='c_cap_OMe'),
            pytest.param('ac-A-_OBn', id='c_cap_OBn'),
            pytest.param('fmoc-G-G-OEt', id='c_cap_OEt_multi'),
        ],
    )
    def test_capped_peptide_partition_roundtrip(self, biln):
        """Amide and ester caps retain product and monomer partitions."""
        self._r(biln)

    # ------------------------------------------------------------------
    # Missing library monomers and explicit synthetic status
    # ------------------------------------------------------------------

    def test_no_question_mark_tokens(self):
        """Registered examples return known monomers without placeholder tokens."""
        peptides = [
            "ac-A-G-V-L-I-P-F-W-M-am",
            "ac-S-T-C-Y-H-D-E-N-Q-am",
            "ac-K-R-am",
            "!1-A-G-V-L-!1",
        ]
        for biln in peptides:
            result, _ = self._r(biln)
            assert "?" not in result, (
                f"Unexpected '?' token in output for {biln!r}: {result!r}"
            )

    def test_library_gap_uses_verified_synthetic_fallback(self, tmp_path, monkeypatch):
        """A library without this residue retains exact chemistry locally."""
        from pyPept.interfaces.cli_monomer import register_monomer
        from pyPept.smiles import convert_smiles

        library = tmp_path / 'alanine-only.sdf'
        register_monomer('N[C@@H](C)C(=O)O', symbol='OnlyAla', sdf_path=library)
        monkeypatch.setenv('CABILN_MONOMER_LIBRARY', str(library))
        smiles = 'NCC(=O)O'
        result = convert_smiles(smiles)
        assert result.recognition_status == 'partial'
        assert len(result.assignments) == 1
        assert not result.assignments[0].recognized
        assert result.details[0][1] == 0
        rebuilt = Molecule(Sequence(result.cabiln)).get_molecule(fmt='ROMol')
        assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

    # ------------------------------------------------------------------
    # Edge-case tests 14-20, 22-30
    # ------------------------------------------------------------------

    def test_explicit_h_amide_smiles(self):
        """Explicit H on backbone amide N/C atoms does not confuse residue matching."""
        from pyPept.smiles import smiles_to_cabiln_core
        smi = 'CC(=O)[NH][C@@H](C)C(=O)[NH2]'  # ac-A-am with explicit H
        result, details = smiles_to_cabiln_core(smi)
        assert 'A' in result, f"Expected A residue in output, got {result!r}"
        assert '?' not in result, f"Unexpected '?' in {result!r}"

    def test_protonated_sidechain_roundtrip(self):
        """The supplied Lys and Arg sidechain structures survive the round trip."""
        result, _ = self._r('ac-K-R-am')
        assert 'K' in result, f"Lys residue lost: {result!r}"
        assert 'R' in result, f"Arg residue lost: {result!r}"
        assert '?' not in result

    def test_isotope_labels_not_in_output(self):
        """smiles_to_cabiln_core must not leak isotope-label dummy atoms into output."""
        import re as _re
        result, _ = self._r('ac-A-G-V-am')
        assert not _re.search(r'\[\d+\*\]', result), (
            f"Isotope-label dummy found in CABILN output: {result!r}"
        )

    @pytest.mark.parametrize("biln", [
        pytest.param("ac-W-am", id="trp"),
        pytest.param("ac-H-am", id="his"),
        pytest.param("ac-Y-am", id="tyr"),
        pytest.param("ac-F-am", id="phe"),
        pytest.param("ac-P-am", id="pro"),
    ])
    def test_aromatic_sidechain_roundtrip(self, biln):
        """Aromatic and cyclic sidechains (W/H/Y/F/P) survive the SMILES roundtrip."""
        expected_abbr = biln[3]  # 'W', 'H', 'Y', 'F', or 'P'
        result, _ = self._r(biln)
        assert expected_abbr in result, (
            f"Expected {expected_abbr!r} in roundtrip output for {biln!r}, got {result!r}"
        )
        assert '?' not in result

    def test_cabiln_to_bracket_idempotent(self):
        """Applying cabiln_to_bracket twice returns the same string (idempotent)."""
        from pyPept.sequence import cabiln_to_bracket
        original = 'ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am'
        once = cabiln_to_bracket(original)
        twice = cabiln_to_bracket(once)
        assert once == twice, (
            f"cabiln_to_bracket not idempotent:\n  first:  {once!r}\n  second: {twice!r}"
        )

    def test_boc_ncap_roundtrip(self):
        """boc N-cap on a peptide with a disulfide crosslink round-trips correctly."""
        result, _ = self._r('boc-C.!1(4,4)-A-G-C.!1(4,4)-am')
        assert result.lower().startswith('boc-'), f"boc cap lost: {result!r}"
        assert '.!1(4,4)' in result, f"Disulfide annotation lost: {result!r}"
        assert '?' not in result

    def test_fmoc_ncap_with_scaffold_roundtrip(self):
        """fmoc N-cap with TATA scaffold annotation survives SMILES roundtrip."""
        result, _ = self._r('fmoc-C.!1(4,4)-A-A-C.!2(4,5)-am%TATA.!1.!2')
        assert result.startswith('fmoc-'), f"fmoc cap lost: {result!r}"
        assert '%TATA' in result, f"TATA scaffold lost: {result!r}"

    def test_n_terminal_proline_assembles(self):
        """N-terminal Pro (secondary amine backbone) assembles and round-trips cleanly."""
        result, _ = self._r('ac-P-A-G-am')
        assert result.startswith('ac-P-'), f"Pro at position 1 lost: {result!r}"
        assert '?' not in result

    def test_nmethyl_aa_chain_roundtrip(self):
        """N-methyl amino acid (meA) in chain round-trips correctly."""
        result, _ = self._r('ac-meA-G-A-am')
        assert 'meA' in result, f"meA residue lost: {result!r}"
        assert '?' not in result

    def test_scaffold_detected_at_most_once(self):
        """smiles_to_cabiln_core annotates a TBMB scaffold exactly once."""
        import re as _re
        result, _ = self._r('ac-C.!1(4,4)-A-C.!2(4,5)-A-C.!3(4,6)-am%TBMB.!1.!2.!3')
        count = len(_re.findall(r'%TBMB', result))
        assert count == 1, f"Expected exactly 1 %TBMB, found {count} in {result!r}"

    def test_scaffold_patterns_1arm_halogen_detected(self):
        """One reacted TBMB arm and both remaining brominated arms are preserved."""
        result, _ = self._r('ac-C.!1(4,4)-G-A-K-C-E-L-F-C-am%TBMB.!1')
        assert '%TBMB' in result, (
            f"1-arm TBMB with Br leaving groups must carry scaffold annotation: {result!r}"
        )
        assert result.split('%TBMB.')[1].count('!') == 1, (
            f"Expected single-arm annotation %%TBMB.!1, got {result!r}"
        )

    def test_scaffold_patterns_1arm_tata_detected(self):
        """One reacted TATA arm and both unreacted arms retain their assembled structure."""
        result, _ = self._r('ac-C.!1(4,4)-A-A-am%TATA.!1')
        assert '%TATA' in result, (
            f"1-arm TATA must carry scaffold annotation: {result!r}"
        )
        assert result.split('%TATA.')[1].count('!') == 1, (
            f"Expected single-arm annotation %TATA.!1, got {result!r}"
        )

    def test_tbmb_mol_not_annotated_as_tata(self):
        """A TBMB-crosslinked peptide is labelled %TBMB, never %TATA."""
        result, _ = self._r('ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2')
        assert '%TBMB' in result, f"Expected %TBMB in {result!r}"
        assert '%TATA' not in result, f"Spurious %TATA annotation in {result!r}"

    # ── Short capped chains ──────────────────────────────────────────────────

    @pytest.mark.parametrize("biln", [
        pytest.param("ac-G-G-am",   id="ac_G_G_am"),
        pytest.param("ac-G-A-am",   id="ac_G_A_am"),
        pytest.param("ac-G-G-G-am", id="ac_G_G_G_am"),
        pytest.param("fmoc-G-am",   id="fmoc_G_am"),
    ])
    def test_multi_amide_backbone_less_chain_roundtrip(self, biln):
        """Short capped peptides retain all monomer occurrences and attachment boundaries."""
        result, _ = self._r(biln)

    def test_large_peptide_performance(self):
        """20-residue peptide round-trips in under 10 seconds."""
        import time as _time
        residues = ['A', 'G', 'V', 'L', 'I', 'P', 'F', 'W', 'M', 'S'] * 2
        biln = 'ac-' + '-'.join(residues) + '-am'
        t0 = _time.time()
        result, _ = self._r(biln)
        elapsed = _time.time() - t0
        assert elapsed < 10.0, f"20-mer round-trip took {elapsed:.1f}s (limit 10s)"
        assert '?' not in result, f"Unexpected '?' in 20-mer output: {result!r}"


class TestGLP1DrugSMILES:
    """Molecular-fidelity regressions for the six named GLP-1 drug fixtures."""

    @staticmethod
    def _convert_smiles(smi: str):
        from pyPept.smiles import convert_smiles
        return convert_smiles(smi)

    @pytest.mark.parametrize("drug,smiles,expected_cabiln,n_res", [
        pytest.param(
            "semaglutide",
            # Wikipedia — C187H291N45O59, MW 4113.6
            (
                "CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCC(=O)O)NC(=O)"
                "[C@H](CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@H](NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)"
                "NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(=O)O)"
                "NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)"
                "NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)"
                "NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)NC(=O)[C@@H](N)Cc1c[nH]cn1)"
                "[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)"
                "C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)"
                "N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)C(C)C"
            ),
            "H-Aib-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K.[AEEA(4,2).AEEA(1,2).E_g(1,2).C18FA(1,2)]-E-F-I-A-W-L-V-R-G-R-G",
            31,
            id="semaglutide",
        ),
        pytest.param(
            "liraglutide",
            # Wikipedia — C172H265N43O51, MW 3751.3
            (
                "CCCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)NCCCC[C@H](NC(=O)[C@H](C)NC(=O)[C@H](C)"
                "NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(C)C)"
                "NC(=O)[C@H](CC1=CC=C(O)C=C1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@@H]"
                "(NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC1=CC=CC=C1)"
                "NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](C)NC(=O)[C@@H](N)"
                "CC1=CN=CN1)[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](CCC(O)=O)C(=O)"
                "N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C)C(=O)"
                "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)"
                "N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(O)=O)C(O)=O"
            ),
            "H-A-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K.[E(4,4).Pal(1,2)]-E-F-I-A-W-L-V-R-G-R-G",
            31,
            id="liraglutide",
        ),
        pytest.param(
            "exenatide",
            # Wikipedia — C184H282N50O60S, MW 4186.6
            (
                "[H]/N=C(\\N)/NCCC[C@@H](C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1ccccc1)C(=O)"
                "N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2c[nH]c3c2cccc3)"
                "C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)"
                "NCC(=O)N4CCC[C@H]4C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)"
                "C(=O)N5CCC[C@H]5C(=O)N6CCC[C@H]6C(=O)N7CCC[C@H]7C(=O)N[C@@H](CO)C(=O)N)"
                "NC(=O)[C@H](C(C)C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)"
                "NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CCCCN)"
                "NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)"
                "NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](Cc8ccccc8)NC(=O)[C@H]([C@@H](C)O)"
                "NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](Cc9cnc[nH]9)N"
            ),
            "H-G-E-G-T-F-T-S-D-L-S-K-Q-M-E-E-E-A-V-R-L-F-I-E-W-L-K-N-G-G-P-S-S-G-A-P-P-P-S-am",
            39,
            id="exenatide",
        ),
        pytest.param(
            "lixisenatide",
            # PubChem CID 90472060 — MW 4858.6
            (
                "CC[C@H](C)[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)"
                "C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)"
                "NCC(=O)N3CCC[C@H]3C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)"
                "C(=O)N4CCC[C@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCCN)"
                "C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)"
                "C(=O)N[C@@H](CCCCN)C(=O)N)NC(=O)[C@H](CC6=CC=CC=C6)NC(=O)[C@H](CC(C)C)"
                "NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](C(C)C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)"
                "NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCC(=O)N)"
                "NC(=O)[C@H](CCCCN)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)"
                "NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC7=CC=CC=C7)NC(=O)"
                "[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](CC8=CNC=N8)N"
            ),
            "H-G-E-G-T-F-T-S-D-L-S-K-Q-M-E-E-E-A-V-R-L-F-I-E-W-L-K-N-G-G-P-S-S-G-A-P-P-S-K-K-K-K-K-K-am",
            44,
            id="lixisenatide",
        ),
        pytest.param(
            "tirzepatide",
            # PubChem CID 166567236 — C225H348N48O68, MW 4813.5
            (
                "CC[C@H](C)[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H]"
                "(CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@@H](C(=O)O)NC(=O)CCCCCCCCCCCCCCCCCCC(=O)O)"
                "C(=O)N[C@@H](C)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](C(C)C)C(=O)"
                "N[C@@H](CCC(=O)N)C(=O)N[C@@H](CC2=CNC3=CC=CC=C32)C(=O)N[C@@H](CC(C)C)"
                "C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C)C(=O)NCC(=O)NCC(=O)N4CCC[C@H]4"
                "C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)C(=O)N5CCC[C@H]5"
                "C(=O)N6CCC[C@H]6C(=O)N7CCC[C@H]7C(=O)N[C@@H](CO)C(=O)N)NC(=O)[C@H](CCCCN)"
                "NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)C(C)(C)NC(=O)[C@H]([C@@H](C)CC)"
                "NC(=O)[C@H](CO)NC(=O)[C@H](CC8=CC=C(C=C8)O)NC(=O)[C@H](CC(=O)O)NC(=O)"
                "[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC9=CC=CC=C9)NC(=O)"
                "[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)"
                "NC(=O)[C@H](CC1=CC=C(C=C1)O)N"
            ),
            "Y-Aib-E-G-T-F-T-S-D-Y-S-I-Aib-L-D-K-I-A-Q-K.[AEEA(4,2).AEEA(1,2).E_g(1,2).C20FA(1,2)]-A-F-V-Q-W-L-I-A-G-G-P-S-S-G-A-P-P-P-S-am",
            39,
            id="tirzepatide",
        ),
        pytest.param(
            "retatrutide",
            # PubChem CID 171390338 — MW 4731.4 (triple GIP/GLP-1/GCG agonist)
            (
                "CC[C@H](C)C(C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CC=C(C=C1)O)C(=O)"
                "N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)NCC(=O)"
                "NCC(=O)N2CCC[C@H]2C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)"
                "C(=O)N3CCC[C@H]3C(=O)N4CCC[C@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](CO)C(=O)N)"
                "NC(=O)[C@H](CC6=CC=CC=C6)NC(=O)[C@H](C)NC(=O)C(C)(C)NC(=O)[C@H](CCC(=O)N)"
                "NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)COCCOCCNC(=O)CC[C@@H](C(=O)O)"
                "NC(=O)CCCCCCCCCCCCCCCCCCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(=O)O)"
                "NC(=O)[C@H](CC(C)C)NC(=O)[C@@](C)(CC(C)C)NC(=O)C([C@@H](C)CC)NC(=O)"
                "[C@H](CO)NC(=O)[C@H](CC7=CC=C(C=C7)O)NC(=O)[C@H](CC(=O)O)NC(=O)"
                "[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC8=CC=CC=C8)NC(=O)"
                "[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)N)NC(=O)C(C)(C)"
                "NC(=O)[C@H](CC9=CC=C(C=C9)O)N"
            ),
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am",
            39,
            id="retatrutide",
        ),
    ])
    def test_glp1_drug_smiles_roundtrip(self, drug, smiles, expected_cabiln, n_res):
        """Preserve the supplied molecular graph and every defined stereo center.

        The historical expected strings identify useful library decompositions,
        but some silently change a histidine tautomer or omit defined stereo.
        Exact isomeric SMILES equality is required for fully specified fixtures.
        Retatrutide omits two alpha configurations: library inference is allowed
        there only if connectivity and every supplied configuration are retained.
        """
        result = self._convert_smiles(smiles)
        cabiln, details = result.cabiln, result.details
        source = Chem.MolFromSmiles(smiles)
        rebuilt = Molecule(Sequence(cabiln)).get_molecule(fmt='ROMol')
        assert all(abbr != '?' for abbr, _, _ in details)
        if drug != 'retatrutide':
            assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(source)
        else:
            source_flat, rebuilt_flat = Chem.Mol(source), Chem.Mol(rebuilt)
            Chem.RemoveStereochemistry(source_flat)
            Chem.RemoveStereochemistry(rebuilt_flat)
            assert Chem.MolToSmiles(rebuilt_flat) == Chem.MolToSmiles(source_flat)
            assert rebuilt.HasSubstructMatch(source, useChirality=True)
        # Molecular identity and editable residue decomposition are independent
        # requirements. A whole-component synthetic fallback fails this check.
        assert len(details) == n_res
        assert result.inferred_stereo == (drug == 'retatrutide')
        if drug in {'semaglutide', 'exenatide', 'lixisenatide'}:
            synthetic = [abbr for abbr, _, _ in details if abbr.startswith('<')]
            assert len(synthetic) == 1
            ordinary_name = 'R' if drug == 'exenatide' else 'H'
            _assert_same_monomer_partition(expected_cabiln, cabiln.replace(synthetic[0], ordinary_name))
            assert result.synthetic_components == (0,)
            assert any('local synthetic' in message for message in result.warnings)
        else:
            assert result.synthetic_components == ()
        if drug in {'liraglutide', 'tirzepatide'}:
            _assert_same_monomer_partition(expected_cabiln, cabiln)
