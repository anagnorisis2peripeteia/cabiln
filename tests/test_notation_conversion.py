"""Notation conversions and molecular equivalence of alternate forms."""

import os
import sys
import warnings

import pytest
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.molecule import Molecule
from pyPept.sequence import (
    Sequence,
    biln_to_cabiln,
)

from _chemistry_oracles import (
    _romol,
    _smiles,
)


class TestRoundTrips:
    """Round-trip tests: BILN, HELM, and FASTA all produce correct molecules."""

    # ------------------------------------------------------------------ #
    # biln_to_cabiln utility
    # ------------------------------------------------------------------ #

    def test_biln_to_cabiln_symmetric(self):
        """Old BILN R3 sidechain crosslink remaps to pyPept slot 4 (backbone_n_mod at slot 3)."""
        result = biln_to_cabiln('C(1,3)-A-A-A-C(1,3)')
        assert '.!1(4,4)' in result   # R3 → slot 4 (thiol for Cys)
        assert result.count('.!1') == 2

    def test_biln_to_cabiln_asymmetric(self):
        """Old BILN asymmetric crosslink: R3 slots remapped to pyPept slot 4."""
        result = biln_to_cabiln('K(1,4)-A-A-A-D(1,3)')
        assert '.!1(4,4)' in result   # K R4 stays 4; D R3 → slot 4 (sidechain carboxyl)
        assert 'D.!1' in result       # second endpoint: implicit

    def test_biln_to_cabiln_no_op(self):
        """biln_to_cabiln leaves strings with no old crosslinks unchanged."""
        s = 'fmoc-A-G-K-am'
        assert biln_to_cabiln(s) == s

    def test_biln_to_cabiln_two_bonds(self):
        """Two independent crosslinks both get converted."""
        result = biln_to_cabiln('C(1,3)-A-K(2,4)-A-C(1,3)-A-D(2,3)')
        assert '.!1' in result
        assert '.!2' in result

    # ------------------------------------------------------------------ #
    # Old BILN → explicit fmt='biln' conversion
    # ------------------------------------------------------------------ #

    def test_old_biln_raises_without_fmt(self):
        """Old BILN notation raises ValueError when fmt='biln' not passed."""
        with pytest.raises(ValueError, match='Old BILN'):
            Sequence('C(1,3)-A-A-A-C(1,3)')

    def test_old_biln_crosslink_assembles_with_fmt(self):
        """Old BILN disulfide assembles correctly when fmt='biln' is passed."""
        seq = Sequence('C(1,3)-A-A-A-C(1,3)', fmt='biln')
        assert seq is not None
        from pyPept.molecule import Molecule
        from rdkit import Chem
        romol = Molecule(seq).get_molecule(fmt='ROMol')
        assert romol is not None
        assert Chem.MolToSmiles(romol) != ''

    # ------------------------------------------------------------------ #
    # HELM round-trip
    # ------------------------------------------------------------------ #

    def test_helm_linear_round_trip(self):
        """HELM linear → Converter.get_biln() → Sequence → molecule works."""
        from pyPept.converter import Converter
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            b = Converter(helm='PEPTIDE1{[ac].D.T.H.F.E.I.A.[am]}$$$$V2.0')
            biln = b.get_biln()
            romol = _romol(biln)
        assert romol is not None

    def test_helm_crosslink_round_trip(self):
        """HELM crosslink → Converter.get_biln() → Sequence → molecule works."""
        from pyPept.converter import Converter
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            b = Converter(helm='PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0')
            biln = b.get_biln()
            assert '.!1' in biln, f"Expected CABILN notation in: {biln}"
            romol = _romol(biln)
        assert romol is not None

    def test_helm_crosslink_biln_contains_cabiln_notation(self):
        """Converter.get_biln() emits CABILN .!n notation, not old (bid,rg)."""
        from pyPept.converter import Converter
        import re
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            b = Converter(helm='PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0')
            biln = b.get_biln()
        # Must not contain bare integer crosslink notation
        assert not re.search(r'(?<!.[\w])[A-Za-z]\w*\(\d+,\d+\)', biln), \
            f"Old BILN notation present in: {biln}"
        assert '.!1' in biln

    # ------------------------------------------------------------------ #
    # FASTA round-trip
    # ------------------------------------------------------------------ #

    def test_fasta_round_trip(self):
        """FASTA single-letter sequence → joined BILN → molecule works."""
        fasta = 'PEPTIDE'
        biln = '-'.join(list(fasta))
        romol = _romol(biln)
        assert romol is not None

    def test_fasta_all_natural_aas(self):
        """All 20 standard amino acids assemble via FASTA-style joined BILN."""
        fasta = 'ACDEFGHIKLMNPQRSTVWY'
        biln = '-'.join(list(fasta))
        romol = _romol(biln)
        assert romol is not None


class TestNotationConversion:
    """Branch <-> bracket notation conversion round-trips.

    Tests cabiln_to_branch() and cabiln_to_bracket() with:
      - Positional branch matching (no !n tags on branches)
      - Continuation monomers carry their R-group pairs
      - Mixed notation (multiple branches, branches + crosslinks)
      - Identity pass-through (no brackets / no %)
      - Multi-monomer branches, single-monomer branches

    String-equality assertions verify the notation form.
    SMILES-equality assertions (via _smiles) verify molecular identity
    for roundtrip tests where both ends assemble.
    """

    @staticmethod
    def _smiles(cabiln):
        """Assemble CABILN and return canonical SMILES.

        Positional % branch notation (no !n crosslink marker on the branch
        monomers) cannot be parsed directly by Sequence — it must be in
        bracket [...] form.  This helper detects that case and converts
        automatically so callers don't need to care.
        """
        from pyPept.sequence import Sequence, cabiln_to_bracket
        from pyPept.molecule import Molecule
        from rdkit.Chem import MolToSmiles
        if '%' in cabiln:
            branch_parts = cabiln.split('%')[1:]
            # Convert unless EVERY branch segment uses an explicit !n crosslink
            # marker — mixed forms (some crosslink, some positional) also need
            # conversion because Sequence can't handle positional .(r,r) syntax.
            all_crosslink = all('!' in p for p in branch_parts)
            if not all_crosslink:
                cabiln = cabiln_to_bracket(cabiln)
        romol = Molecule(Sequence(cabiln)).get_molecule(fmt="ROMol")
        return MolToSmiles(romol)

    # ------------------------------------------------------------------
    # Branch -> Bracket -> Branch  (% notation is canonical input)
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("label,branch,expected_bracket", [
        ("lipid_linker",
         "ac-K.!1(4,4)-G-am%C20FA-AEEA(2,1)-E_g(2,1).!1",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"),

        ("NH2_branch",
         "ac-A-D.(4,1)-G-am%G-A(2,1)-am(2,1)",
         "ac-A-D.[G(4,1).A(2,1).am(2,1)]-G-am"),

        ("COOH_branch",
         "A-D.(4,1)-A-am%G-A(2,1)-am(2,1)",
         "A-D.[G(4,1).A(2,1).am(2,1)]-A-am"),

        ("disulfide_branch",
         "ac-C.(4,4)-A-G-am%C-A(2,1)-am(2,1)",
         "ac-C.[C(4,4).A(2,1).am(2,1)]-A-G-am"),

        ("single_monomer",
         "ac-K.(4,2)-G-am%A",
         "ac-K.[A(4,2)]-G-am"),

        ("unannotated_defaults_2_1",
         "ac-D.(4,1)-G-am%G-A-am",
         "ac-D.[G(4,1).A(2,1).am(2,1)]-G-am"),
    ])
    def test_branch_to_bracket(self, label, branch, expected_bracket):
        from pyPept.sequence import cabiln_to_bracket
        result = cabiln_to_bracket(branch)
        assert result == expected_bracket, (
            f"{label}: expected {expected_bracket!r}, got {result!r}"
        )

    @pytest.mark.parametrize("label,branch", [
        ("NH2_branch",
         "ac-A-D.!1(4,1)-G-am%G.!1-A-am"),
        ("COOH_branch",
         "A-D.!1(4,1)-A-am%G.!1-A-am"),
        ("disulfide_branch",
         "ac-C.!1(4,4)-A-G-am%C.!1-A-am"),
    ])
    def test_branch_roundtrip(self, label, branch):
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = cabiln_to_bracket(branch)
        back = cabiln_to_branch(bracket)
        assert back == branch, (
            f"{label}: roundtrip failed\n"
            f"  branch  -> bracket: {bracket!r}\n"
            f"  bracket -> branch:  {back!r}"
        )

    def test_crosslink_branch_to_bracket_roundtrip(self):
        """Crosslink !n branches convert to (1,2) bracket and back."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-K.!1(4,4)-G-am%C20FA-AEEA-E_g.!1"
        bracket = cabiln_to_bracket(branch)
        assert bracket == "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        back = cabiln_to_branch(bracket)
        assert back == "ac-K.!1(4,4)-G-am%C20FA-AEEA-E_g.!1"
        assert self._smiles(bracket) == self._smiles(back)

    # ------------------------------------------------------------------
    # Bracket -> Branch -> Bracket  (bracket notation is canonical input)
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("label,bracket,expected_branch", [
        ("bracket_COOH",
         "A-D.[G(4,1).A(2,1).am(2,1)]-A-am",
         "A-D.!1(4,1)-A-am%G.!1-A-am"),

        ("bracket_NH2",
         "ac-A-D.[G(4,1).A(2,1).am(2,1)]-G-am",
         "ac-A-D.!1(4,1)-G-am%G.!1-A-am"),

        ("bracket_disulfide",
         "ac-C.[C(4,4).A(2,1).am(2,1)]-A-G-am",
         "ac-C.!1(4,4)-A-G-am%C.!1-A-am"),

        ("bracket_single_monomer",
         "ac-K.[A(4,2)]-G-am",
         "ac-K.[A(4,2)]-G-am"),
    ])
    def test_bracket_to_branch(self, label, bracket, expected_branch):
        from pyPept.sequence import cabiln_to_branch
        result = cabiln_to_branch(bracket)
        assert result == expected_branch, (
            f"{label}: expected {expected_branch!r}, got {result!r}"
        )

    @pytest.mark.parametrize("label,bracket", [
        ("bracket_COOH",
         "A-D.[G(4,1).A(2,1).am(2,1)]-A-am"),
        ("bracket_NH2",
         "ac-A-D.[G(4,1).A(2,1).am(2,1)]-G-am"),
        ("bracket_disulfide",
         "ac-C.[C(4,4).A(2,1).am(2,1)]-A-G-am"),
    ])
    def test_bracket_roundtrip(self, label, bracket):
        from pyPept.sequence import cabiln_to_branch, cabiln_to_bracket
        branch = cabiln_to_branch(bracket)
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"{label}: roundtrip failed\n"
            f"  bracket -> branch:  {branch!r}\n"
            f"  branch  -> bracket: {back!r}"
        )

    # ------------------------------------------------------------------
    # Identity / pass-through
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("notation", [
        "ac-A-G-K-am",
        "A-G",
        "am",
        "",
    ])
    def test_no_branch_passthrough(self, notation):
        from pyPept.sequence import cabiln_to_branch, cabiln_to_bracket
        assert cabiln_to_branch(notation) == notation
        assert cabiln_to_bracket(notation) == notation

    # ------------------------------------------------------------------
    # Mixed notation — multiple branches (positional matching)
    # ------------------------------------------------------------------

    def test_two_branches_roundtrip(self):
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-C.!1(4,4)-A-K.!2(4,1)-G-am%C.!1-A-am%G.!2-am"
        bracket = cabiln_to_bracket(branch)
        assert "[" in bracket and "]" in bracket
        back = cabiln_to_branch(bracket)
        assert back == branch, (
            f"Two-branch roundtrip failed:\n"
            f"  {branch!r} -> {bracket!r} -> {back!r}"
        )

    def test_two_brackets_roundtrip(self):
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = (
            "ac-C.[C(4,4).A(2,1).am(2,1)]"
            "-A-D.[G(4,1).am(2,1)]-G-am"
        )
        branch = cabiln_to_branch(bracket)
        assert "%" in branch, f"Should have branch separators: {branch}"
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"Two-bracket roundtrip failed:\n"
            f"  {bracket!r} -> {branch!r} -> {back!r}"
        )

    # ------------------------------------------------------------------
    # Crosslink + branch coexistence
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Unannotated continuation defaults to (2,1)
    # ------------------------------------------------------------------

    def test_unannotated_defaults_2_1(self):
        from pyPept.sequence import cabiln_to_bracket
        branch = "A-D.(4,1)-am%G-A"
        bracket = cabiln_to_bracket(branch)
        assert "(2,1)" in bracket, (
            f"Unannotated continuation should default to (2,1): {bracket}"
        )

    # ------------------------------------------------------------------
    # (1,2) brackets convert to reversed crosslink branch
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("label, bracket, expected_branch", [
        ("lipid_linker_12",
         "K.[C20FA(4,4).AEEA(1,2).E_g(1,2)]-G-am",
         "K.!1(4,4)-G-am%E_g-AEEA-C20FA.!1"),
        ("two_monomer_12",
         "ac-A.[C(4,4).G(1,2)]-am",
         "ac-A.!1(4,4)-am%G-C.!1"),
        ("lipid_correct",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am",
         "ac-K.!1(4,4)-G-am%C20FA-AEEA-E_g.!1"),
    ])
    def test_12_bracket_converts_to_reversed_branch(self, label, bracket, expected_branch):
        from pyPept.sequence import cabiln_to_branch
        result = cabiln_to_branch(bracket)
        assert result == expected_branch, (
            f"{label}: expected {expected_branch!r}, got {result!r}"
        )

    @pytest.mark.parametrize("label, branch, expected_bracket", [
        (
            "crosslink_lipid_to_bracket",
            "K.!1(4,4)-G-am%E_g-AEEA(2,1)-C20FA(2,1).!1",
            "K.[C20FA(4,4).AEEA(1,2).E_g(1,2)]-G-am",
        ),
        (
            "crosslink_two_monomer_to_bracket",
            "ac-A.!1(4,4)-am%G-C(2,1).!1",
            "ac-A.[C(4,4).G(1,2)]-am",
        ),
    ])
    def test_crosslink_branch_to_bracket(self, label, branch, expected_bracket):
        """Crosslink !n branches still convert to bracket (one-way)."""
        from pyPept.sequence import cabiln_to_bracket
        result = cabiln_to_bracket(branch)
        assert result == expected_bracket, (
            f"{label}: branch->bracket mismatch\n"
            f"  input:    {branch!r}\n"
            f"  expected: {expected_bracket!r}\n"
            f"  got:      {result!r}"
        )

    # ------------------------------------------------------------------
    # Midpoint anchor — anchor in the middle of the branch chain
    # Bracket: K.[E_g(4,4).AEEA(2,1).C(2,1).A(1,2).G(1,2)]-am
    #   E_g is anchor, AEEA-C grow C-terminal (2,1), A-G grow N-terminal (1,2)
    # Branch: anchor at midpoint → need crosslink for the N-terminal arm
    # ------------------------------------------------------------------

    def test_12_lipid_bracket_converts(self):
        """(1,2) lipid linker bracket converts to reversed crosslink branch."""
        from pyPept.sequence import cabiln_to_branch
        bracket = "K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        result = cabiln_to_branch(bracket)
        assert result == "K.!1(4,4)-G-am%C20FA-AEEA-E_g.!1", f"Got: {result}"

    # ------------------------------------------------------------------
    # Complex round-trips
    # ------------------------------------------------------------------

    def test_cyclic_plus_lipid_branch_roundtrip(self):
        """Head-to-tail cyclic peptide with lipid branch: both notation directions."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "!1-A-K.!2(4,4)-G-A-!1%C20FA-AEEA-E_g.!2"
        expected_bracket = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        bracket = cabiln_to_bracket(branch)
        assert bracket == expected_bracket, f"branch->bracket: {bracket!r}"
        back = cabiln_to_branch(bracket)
        assert back == branch, (
            f"bracket->branch roundtrip failed:\n"
            f"  {branch!r} -> {bracket!r} -> {back!r}"
        )
        assert self._smiles(bracket) == self._smiles(back)


    def test_long_chain_two_positional_branches_roundtrip(self):
        """9-residue chain with two positional branches — D (isopeptide) and E (isopeptide)."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-A-G-D.!1(4,1)-A-A-E.!2(4,1)-G-am%G.!1-A-am%A.!2-G-am"
        expected_bracket = (
            "ac-A-G-D.[G(4,1).A(2,1).am(2,1)]"
            "-A-A-E.[A(4,1).G(2,1).am(2,1)]-G-am"
        )
        bracket = cabiln_to_bracket(branch)
        assert bracket == expected_bracket, f"branch->bracket: {bracket!r}"
        back = cabiln_to_branch(bracket)
        assert back == branch, (
            f"bracket->branch roundtrip failed:\n"
            f"  {branch!r} -> {bracket!r} -> {back!r}"
        )

    def test_mixed_bracket_types_roundtrip(self):
        """Single-monomer bracket (passthrough) + lipid (1,2) bracket coexist."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        # [A(4,1)] single-monomer stays as-is; lipid converts to crosslink branch
        bracket = "ac-K.[A(4,1)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am"
        expected_branch = "ac-K.[A(4,1)]-G-K.!1(4,4)-am%C20FA-AEEA-E_g.!1"
        branch = cabiln_to_branch(bracket)
        assert branch == expected_branch, f"bracket->branch: {branch!r}"
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"branch->bracket roundtrip failed:\n"
            f"  {bracket!r} -> {branch!r} -> {back!r}"
        )

    # ------------------------------------------------------------------
    # Comprehensive parametrized round-trip suite
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("label,bracket,expected_branch", [
        ("two_lipid_brackets",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am",
         "ac-K.!1(4,4)-G-K.!2(4,4)-am%C20FA-AEEA-E_g.!1%C20FA-AEEA-E_g.!2"),
        ("two_aeea_peg",
         "ac-K.[AEEA(4,2).AEEA(1,2)]-G-am",
         "ac-K.!1(4,2)-G-am%AEEA-AEEA-!1"),
        ("three_aeea_peg",
         "ac-K.[AEEA(4,2).AEEA(1,2).AEEA(1,2)]-G-am",
         "ac-K.!1(4,2)-G-am%AEEA-AEEA-AEEA-!1"),
        ("four_residue_branch",
         "ac-D.[G(4,1).A(2,1).G(2,1).am(2,1)]-G-am",
         "ac-D.!1(4,1)-G-am%G.!1-A-G-am"),
        ("no_caps",
         "D.[G(4,1).A(2,1).am(2,1)]-G-A",
         "D.!1(4,1)-G-A%G.!1-A-am"),
        ("cyclic_plus_plain",
         "!1-A-D.[G(4,1).am(2,1)]-G-A-!1",
         "!1-A-D.!2(4,1)-G-A-!1%G.!2-am"),
        ("semaglutide_like",
         "His-Aib-Glu-Gly-Thr-Phe-Thr-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am",
         "His-Aib-Glu-Gly-Thr-Phe-Thr-K.!1(4,4)-am%C20FA-AEEA-E_g.!1"),
        ("c_terminal_marker_normalised",
         "ac-K.[G(4,2).G(1,2)]-G-am",
         "ac-K.!1(4,2)-G-am%G-G-!1"),
        ("n_terminal_marker_normalised",
         "ac-D.[G(4,1).G(2,1).am(2,1)]-G-am",
         "ac-D.!1(4,1)-G-am%G.!1-G-am"),
        # Exotic / mixed-type cases
        ("mixed_lipid_isopeptide",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-D.[G(4,1).A(2,1).am(2,1)]-am",
         "ac-K.!1(4,4)-A-D.!2(4,1)-am%C20FA-AEEA-E_g.!1%G.!2-A-am"),
        ("long_isopeptide_arm",
         "ac-E.[G(4,1).G(2,1).A(2,1).A(2,1).am(2,1)]-G-am",
         "ac-E.!1(4,1)-G-am%G.!1-G-A-A-am"),
        ("cyclic_plus_two_residue_arm",
         "!1-A-D.[G(4,1).A(2,1).am(2,1)]-G-A-!1",
         "!1-A-D.!2(4,1)-G-A-!1%G.!2-A-am"),
        ("two_isopeptide_one_lipid",
         "ac-D.[G(4,1).am(2,1)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-D.[A(4,1).am(2,1)]-am",
         "ac-D.!1(4,1)-G-K.!2(4,4)-A-D.!3(4,1)-am%G.!1-am%C20FA-AEEA-E_g.!2%A.!3-am"),
        ("ameleu_slot3_arm",
         "ac-D.[aMeLeu(4,3).G(2,1).am(2,1)]-G-am",
         "ac-D.!1(4,3)-G-am%aMeLeu.!1-G-am"),
        ("dual_e_scaffold",
         "ac-E.[G(4,1).A(2,1).am(2,1)]-G-E.[A(4,1).G(2,1).am(2,1)]-am",
         "ac-E.!1(4,1)-G-E.!2(4,1)-am%G.!1-A-am%A.!2-G-am"),
        ("rgd_pendant_arm",
         "ac-G-D.[R(4,1).G(2,1).D(2,1).am(2,1)]-S-am",
         "ac-G-D.!1(4,1)-S-am%R.!1-G-D-am"),
        ("cyclic_isopeptide_lipid",
         "!1-A-D.[G(4,1).am(2,1)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-!1",
         "!1-A-D.!2(4,1)-G-K.!3(4,4)-A-!1%G.!2-am%C20FA-AEEA-E_g.!3"),
        ('three_positional_branches',
         'ac-D.[G(4,1).am(2,1)]-A-D.[A(4,1).am(2,1)]-A-D.[G(4,1).G(2,1).am(2,1)]-am',
         'ac-D.!1(4,1)-A-D.!2(4,1)-A-D.!3(4,1)-am%G.!1-am%A.!2-am%G.!3-G-am'),
    ])
    def test_bracket_to_branch_comprehensive(self, label, bracket, expected_branch):
        """bracket -> branch: exact string + roundtrip + independent SMILES equivalence."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        result = cabiln_to_branch(bracket)
        assert result == expected_branch, (
            f"{label}: bracket->branch\n"
            f"  expected: {expected_branch!r}\n"
            f"  got:      {result!r}"
        )
        back = cabiln_to_bracket(result)
        assert back == bracket, (
            f"{label}: branch->bracket roundtrip\n"
            f"  expected: {bracket!r}\n"
            f"  got:      {back!r}"
        )
        assert self._smiles(bracket) == self._smiles(result), (
            f"{label}: bracket and branch forms give different molecules"
        )


    def test_retatrutide_bracket_roundtrip(self):
        """Retatrutide 39AA in bracket form: bracket -> branch -> bracket."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
            "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
        )
        expected_branch = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.!1(4,2)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%C20FA-E_g-AEEA-!1"
        )
        branch = cabiln_to_branch(bracket)
        assert branch == expected_branch, (
            f"bracket->branch:\n  expected: {expected_branch!r}\n  got: {branch!r}"
        )
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"branch->bracket roundtrip:\n  expected: {bracket!r}\n  got: {back!r}"
        )

    # ------------------------------------------------------------------
    # Terminal !n marker forms (hyphen-separated) now handled
    # ------------------------------------------------------------------

    def test_n_terminal_marker_branch_to_bracket(self):
        """N-terminal !1 marker: %!1-G-G-am converts to bracket; normalises to positional on back-trip."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-D.!1(4,1)-G-am%!1-G-G-am"
        bracket = cabiln_to_bracket(branch)
        assert bracket == "ac-D.[G(4,1).G(2,1).am(2,1)]-G-am", (
            f"bracket mismatch: {bracket!r}"
        )
        back = cabiln_to_branch(bracket)
        assert back == "ac-D.!1(4,1)-G-am%G.!1-G-am", f"back form: {back!r}"
        assert self._smiles(bracket) == self._smiles(back), "N-terminal marker: bracket vs normalised branch differ"



    def test_nonstd_slot_e14_bracket_to_branch(self):
        """E(1,4) non-standard γ-Glu slot: bracket → multi-crosslink chain branch."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = "K.[AEEA(4,2).E(1,4).C20FA(1,2)]-am"
        expected_branch = "K.!1(4,2)-am%AEEA.!1.!2(1,4)%E.!2.!3(1,2)%C20FA.!3"
        branch = cabiln_to_branch(bracket)
        assert branch == expected_branch, f"bracket->branch: {branch!r}"
        back = cabiln_to_bracket(branch)
        assert back == bracket, f"branch->bracket roundtrip: {back!r}"

    def test_nonstd_slot_retatrutide_e14_roundtrip(self):
        """Full Retatrutide with E(1,4) bracket: bracket → branch → bracket roundtrip."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.[AEEA(4,2).E(1,4).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
            "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
        )
        expected_branch = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.!1(4,2)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%AEEA.!1.!2(1,4)%E.!2.!3(1,2)%C20FA.!3"
        )
        branch = cabiln_to_branch(bracket)
        assert branch == expected_branch, (
            f"bracket->branch:\n  expected: {expected_branch!r}\n  got: {branch!r}"
        )
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"branch->bracket roundtrip:\n  expected: {bracket!r}\n  got: {back!r}"
        )

    def test_nonstd_slot_branch_assembly_matches_bracket(self):
        """E(1,4) branch form assembles to identical SMILES as bracket form."""
        from pyPept.sequence import cabiln_to_branch
        bracket = "K.[AEEA(4,2).E(1,4).C20FA(1,2)]-am"
        branch = cabiln_to_branch(bracket)
        smi_bracket = _smiles(bracket)
        smi_branch = _smiles(branch)
        assert smi_bracket == smi_branch, (
            f"E(1,4) branch/bracket SMILES mismatch\n"
            f"  bracket: {smi_bracket}\n"
            f"  branch:  {smi_branch}"
        )

    @pytest.mark.parametrize("label,bracket,branch", [
        ("disulfide_branch",
         "ac-C.[C(4,4).A(2,1).am(2,1)]-A-G-am",
         "ac-C.(4,4)-A-G-am%C-A-am"),
        ("NH2_branch",
         "ac-A-K.[G(4,2).A(1,2).ac(1,2)]-G-am",
         "ac-A-K.!1(4,2)-G-am%ac-A-G-!1"),
        ("COOH_branch",
         "A-D.[G(4,1).A(2,1).am(2,1)]-A-am",
         "A-D.(4,1)-A-am%G-A-am"),
    ])
    def test_positional_branch_converts_to_assemblable_bracket(
        self, label, bracket, branch
    ):
        """Branch without annotations -> convert to bracket -> assemble, compare."""
        from pyPept.sequence import cabiln_to_bracket
        converted = cabiln_to_bracket(branch)
        assert converted == bracket, f"{label}: conversion mismatch: {converted} != {bracket}"
        smi = self._smiles(bracket)
        assert smi, f"{label}: bracket form should assemble"

    def test_retatrutide_bracket_assembles(self):
        """Full retatrutide with (1,2) lipid linker bracket assembles."""
        cabiln = ("Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
                  "K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
                  "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am")
        smi = self._smiles(cabiln)
        assert smi and len(smi) > 100, f"Retatrutide should produce a large SMILES: {smi}"
