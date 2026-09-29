"""Legacy CABILN lowering and Sequence parsing and validation contracts."""

import os
import sys
import warnings

import pytest
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.sequence import (
    Sequence,
    ValidationReport,
    _preprocess_cabiln,
    _expand_inline_caps,
    colorize_cabiln,
)

from _chemistry_oracles import _romol


class TestInlineAttachments:
    """Tests for .Token(host_r,cap_r) and .!n(r,r) inline notation."""

    # --- _expand_inline_caps unit tests ---

    def test_expand_named_cap_produces_pendant_chain(self):
        """Named cap .ac(2,1) injects crosslink annotation and appends pendant."""
        from pyPept.sequence import _expand_inline_caps
        result, _ = _expand_inline_caps('A-G.ac(2,1)-A')
        # pendant chain appended after a '.'
        assert '.ac(' in result
        # auto bond ID >= 100 injected on host residue
        import re
        ids = re.findall(r'\((\d+),2\)', result)
        assert any(int(i) >= 100 for i in ids), f"No auto ID >= 100 in: {result}"

    def test_expand_crosslink_marker_no_pendant(self):
        """Branch marker .!1(3,3) converts to (!1,3) annotation; no pendant appended."""
        from pyPept.sequence import _expand_inline_caps
        biln = 'A-G.!1(3,3)-A-D.!1(3,3)'
        result, _ = _expand_inline_caps(biln)
        assert '(!1,3)' in result, f"Expected (!1,3) in: {result}"
        # No pendant chain should be appended for !n markers
        parts = result.split('.')
        assert all('!1' not in p.split('(')[0] for p in parts[1:]), \
            f"Unexpected pendant chain with !1 in: {result}"

    def test_expand_chain_separator_unchanged(self):
        """Bare .Token (no parens) remains a chain separator, not an inline cap."""
        from pyPept.sequence import _expand_inline_caps
        result, _ = _expand_inline_caps('A.G')
        assert result == 'A.G', f"Chain separator was altered: {result}"

    def test_expand_multiple_inline_caps(self):
        """Two inline caps on the same residue get distinct auto bond IDs."""
        from pyPept.sequence import _expand_inline_caps
        import re
        result, _ = _expand_inline_caps('Lys.ac(3,1).am(4,1)-G')
        ids = re.findall(r'\((\d+),\d+\)', result)
        numeric_ids = [int(i) for i in ids if not i.startswith('!')]
        # Should have two distinct auto IDs >= 100
        assert len(set(numeric_ids)) >= 2, f"Expected 2 distinct IDs in: {result}"

    def test_expand_inverse_annotation_accepted(self):
        """K.!1(3,1) + G.!1(1,3) — inverse pair is valid, no error."""
        from pyPept.sequence import _expand_inline_caps
        result, br = _expand_inline_caps('K.!1(3,1)-A-G.!1(1,3)-am')
        assert '(!1,3)' in result
        assert '(!1,1)' in result
        assert br['!1'] == 1  # partner rgroup stored from first occurrence

    def test_expand_symmetric_crosslink_accepted(self):
        """K.!1(3,3) + K.!1(3,3) — symmetric crosslink is its own inverse, no error."""
        from pyPept.sequence import _expand_inline_caps
        result, br = _expand_inline_caps('K.!1(3,3)-A-K.!1(3,3)-am')
        assert result.count('(!1,3)') == 2
        assert br['!1'] == 3

    def test_expand_rgroup_conflict_raises(self):
        """K.!1(3,1) + G.!1(2,3) — first says partner uses R1, second declares R2 on self.
        Not the inverse → ValueError."""
        from pyPept.sequence import _expand_inline_caps
        with pytest.raises(ValueError, match='conflict'):
            _expand_inline_caps('K.!1(3,1)-A-G.!1(2,3)-am')

    def test_expand_branch_rgroup_returned(self):
        """branch_rgroup dict maps each !x bond to the partner rgroup from first endpoint."""
        from pyPept.sequence import _expand_inline_caps
        _, br = _expand_inline_caps('C.!1(3,3)-A-C.!1(3,3)-am')
        assert br == {'!1': 3}

    # --- End-to-end assembly tests ---



    def test_crosslink_odd_endpoints_raises(self):
        """Single .!1(3,3) with no matching partner raises ValueError (odd bonds)."""
        with pytest.raises((ValueError, SystemExit)):
            Sequence('K.!1(3,3)-A-am')

    def test_old_biln_raises_without_fmt(self):
        """Old BILN bare-integer crosslink notation raises ValueError unless fmt='biln'."""
        with pytest.raises(ValueError, match='Old BILN'):
            Sequence('K(1,3)-A-K(1,3)-am')

    def test_old_biln_converts_with_fmt_biln(self):
        """Old BILN converts silently when fmt='biln' is passed."""
        seq = Sequence('K(1,3)-A-K(1,3)-am', fmt='biln')
        assert seq is not None

    def test_old_biln_error_mentions_fmt(self):
        """ValueError for old BILN tells user to pass fmt='biln'."""
        with pytest.raises(ValueError, match="fmt='biln'"):
            Sequence('K(1,3)-A-K(1,3)-am')


class TestCABILNPhase2Parser:
    """Unit tests for Phase 2 CABILN notation: %, newline, no-parens .!n, terminals."""

    # --- _preprocess_cabiln ---

    def test_preprocess_strips_whitespace(self):
        assert _preprocess_cabiln('  ac-A-am  ') == 'ac-A-am'

    def test_preprocess_newline_to_percent(self):
        assert _preprocess_cabiln('ac-A-am\nG-G') == 'ac-A-am%G-G'

    def test_preprocess_percent_unchanged(self):
        assert _preprocess_cabiln('ac-A-am%G-G') == 'ac-A-am%G-G'

    def test_preprocess_multiple_newlines_collapsed(self):
        assert _preprocess_cabiln('ac-A-am\n\nG-G') == 'ac-A-am%G-G'

    def test_preprocess_newline_with_spaces(self):
        assert _preprocess_cabiln('ac-A-am \n G-G') == 'ac-A-am%G-G'

    def test_preprocess_strips_leading_trailing_percent(self):
        assert _preprocess_cabiln('%ac-A-am%') == 'ac-A-am'

    # --- _expand_inline_caps: backward-compat (no % separator) ---

    def test_expand_simple_no_change(self):
        result, brg = _expand_inline_caps('ac-A-G-am')
        assert result == 'ac-A-G-am'
        assert brg == {}

    def test_expand_named_cap_unchanged_from_phase1(self):
        result, brg = _expand_inline_caps('fmoc-C.trt(4,2)-am')
        assert 'trt(' in result
        assert '(100,4)' in result
        assert brg == {}

    def test_expand_disulfide_unchanged_from_phase1(self):
        result, brg = _expand_inline_caps('ac-C.!1(4,4)-G-C.!1(4,4)-am')
        assert '(!1,4)' in result
        assert result.count('(!1,4)') == 2
        assert brg == {'!1': 4}

    # --- second occurrence without parens ---

    def test_second_occurrence_no_parens_infers_inverse(self):
        """K.!1(4,2)-G-G.!1 → second G connects via R2 (inverse of R4/R2 → R2/R4)."""
        result, brg = _expand_inline_caps('K.!1(4,2)-G-G.!1')
        assert '(!1,4)' in result   # first occurrence: K uses R4
        assert '(!1,2)' in result   # second occurrence (inferred): G uses R2
        assert brg == {'!1': 2}

    def test_second_occurrence_explicit_inverse_also_works(self):
        """Explicit inverse .!1(2,4) is equivalent to implicit .!1 for second endpoint."""
        result_impl, _ = _expand_inline_caps('K.!1(4,2)-G-G.!1')
        result_expl, _ = _expand_inline_caps('K.!1(4,2)-G-G.!1(2,4)')
        assert result_impl == result_expl

    def test_symmetric_crosslink_no_parens_second(self):
        """Symmetric crosslink .!1(4,4): second .!1 with no parens → also R4."""
        result, brg = _expand_inline_caps('ac-C.!1(4,4)-G-C.!1')
        assert result.count('(!1,4)') == 2
        assert brg == {'!1': 4}

    def test_no_parens_without_prior_occurrence_raises(self):
        """Using .!1 with no parens and no prior .!1(y,z) raises ValueError."""
        with pytest.raises(ValueError, match='first-occurrence'):
            _expand_inline_caps('ac-G.!1-G-am')

    # --- R-group conflict ---

    def test_rgroup_conflict_raises(self):
        """Second .!1(y,z) that is not the inverse of the first raises ValueError."""
        with pytest.raises(ValueError, match='R-group conflict'):
            _expand_inline_caps('ac-C.!1(4,2)-G-C.!1(4,3)-am')

    # --- %-separated branch segments ---

    def test_percent_separator_produces_dot_separated_chains(self):
        """% splits into segments; _expand_inline_caps joins them with '.'."""
        result, _ = _expand_inline_caps('K.!1(4,2)-K-am%G-G.!1')
        assert '.' in result
        chains = result.split('.')
        assert len(chains) == 2

    def test_newline_separator_same_as_percent(self):
        """Newline and % are interchangeable segment separators."""
        res_pct, brg_pct = _expand_inline_caps('K.!1(4,2)-K-am%G-G.!1')
        res_nl,  brg_nl  = _expand_inline_caps('K.!1(4,2)-K-am\nG-G.!1')
        assert res_pct == res_nl
        assert brg_pct == brg_nl

    def test_branch_segment_isopeptide(self):
        """K.!1(4,2)-K-am%G-G-G.!1 → K uses R4; last G uses R2 (isopeptide)."""
        result, brg = _expand_inline_caps('K.!1(4,2)-K-am%G-G-G.!1')
        assert '(!1,4)' in result
        assert '(!1,2)' in result
        assert brg == {'!1': 2}

    # --- terminal markers ---

    def test_c_terminal_marker(self):
        """A-B-C-!1 (C-terminal marker) → last residue annotated with (!1,2)."""
        result, brg = _expand_inline_caps('K.!1(4,2)-K-am%G-G-G-!1')
        assert '(!1,2)' in result
        assert 'G(!1,2)' in result
        assert brg == {'!1': 2}

    def test_n_terminal_marker(self):
        """!1-A-B-C (N-terminal marker) → first residue annotated with (!1,1)."""
        result, brg = _expand_inline_caps('K.!1(4,1)-K-am%!1-G-G-G-am')
        assert '(!1,1)' in result
        assert 'G(!1,1)' in result
        assert brg == {'!1': 1}

    def test_c_terminal_marker_rgroup_mismatch_raises(self):
        """-!1 implies R2 but .!1(4,1) declared partner R1 — should raise."""
        with pytest.raises(ValueError, match='R2|rgroup|attach'):
            _expand_inline_caps('K.!1(4,1)-K-am%G-G-G-!1')

    def test_n_terminal_marker_rgroup_mismatch_raises(self):
        """!1- implies R1 but .!1(4,2) declared partner R2 — should raise."""
        with pytest.raises(ValueError, match='R1|rgroup|attach'):
            _expand_inline_caps('K.!1(4,2)-K-am%!1-G-G-G-am')

    def test_terminal_marker_without_partner_raises(self):
        """A lone terminal marker with no paired endpoint raises on endpoint count."""
        with pytest.raises(ValueError, match='endpoint'):
            _expand_inline_caps('K-K-am%G-G-G-!1')

    def test_three_endpoints_raises(self):
        """Three occurrences of the same bond ID raise ValueError."""
        with pytest.raises(ValueError, match='3rd endpoint|2 endpoint'):
            _expand_inline_caps('C.!1(4,4)-G-C.!1(4,4)-G-C.!1(4,4)-am')

    def test_single_endpoint_raises(self):
        """A bond with only one endpoint raises ValueError."""
        with pytest.raises(ValueError, match='1 endpoint|exactly 2'):
            _expand_inline_caps('ac-C.!1(4,4)-G-am')

    # --- end-to-end Sequence assembly with % notation ---

    def test_sequence_with_percent_branch(self):
        """Sequence accepts %-separated CABILN with a branch pendant chain."""
        seq = Sequence('ac-K.!1(4,2)-G-am%G-G.!1')
        assert seq is not None
        assert seq.s_nmonomers == 6  # ac, K, G, am  +  G, G

    def test_sequence_newline_branch(self):
        """Sequence accepts newline-separated CABILN."""
        seq = Sequence('ac-K.!1(4,2)-G-am\nG-G.!1')
        assert seq is not None
        assert seq.s_nmonomers == 6

    def test_sequence_n_terminal_branch(self):
        """Sequence handles N-terminal branch marker !1-A-B-C."""
        seq = Sequence('ac-K.!1(4,1)-G-am%!1-G-G-am')
        assert seq is not None
        assert seq.s_nmonomers == 7  # ac, K, G, am  +  G, G, am

    def test_sequence_c_terminal_branch(self):
        """Sequence handles C-terminal branch marker A-B-C-!1."""
        seq = Sequence('ac-K.!1(4,2)-G-am%G-G-!1')
        assert seq is not None
        assert seq.s_nmonomers == 6  # ac, K, G, am  +  G, G


class TestSequenceValidate:
    """Sequence.validate() returns a ValidationReport without raising."""

    def test_valid_tripeptide_ok(self):
        report = Sequence.validate('A-G-A')
        assert report.ok is True
        assert report.errors == []
        assert isinstance(report.bonds, list)

    def test_valid_has_backbone_bonds(self):
        """Tripeptide has 2 backbone bonds (A→G and G→A)."""
        report = Sequence.validate('A-G-A')
        assert len(report.bonds) == 2


    def test_valid_build_returns_sequence(self):
        report = Sequence.validate('A-G-A')
        seq = report.build()
        assert isinstance(seq, Sequence)
        assert seq.length() == 3

    def test_invalid_token_not_ok(self):
        """Unknown monomer token → report.ok is False, errors populated."""
        report = Sequence.validate('A-NOTAMONOMER-G')
        assert report.ok is False
        assert len(report.errors) > 0

    def test_invalid_build_raises(self):
        """build() on a failed report raises ValueError."""
        report = Sequence.validate('A-NOTAMONOMER-G')
        with pytest.raises(ValueError, match="Cannot build"):
            report.build()

    def test_valid_capped_peptide(self):
        report = Sequence.validate('ac-A-G-am')
        assert report.ok is True
        assert len(report.bonds) == 3  # ac→A, A→G, G→am

    def test_report_is_validation_report(self):
        report = Sequence.validate('A')
        assert isinstance(report, ValidationReport)


    def test_cabiln_cyclic_peptide_valid(self):
        """Cyclic CABILN !1-A-A-A-A-!1 validates correctly."""
        report = Sequence.validate('!1-A-A-A-A-!1')
        assert report.ok is True
        # 3 backbone bonds + 1 head-to-tail crosslink = 4 bonds
        assert len(report.bonds) == 4

    def test_old_biln_validate_fails_without_fmt(self):
        """Old BILN K(1,3) without fmt='biln' → validate returns not-ok with error."""
        report = Sequence.validate('K(1,3)-A-K(1,3)')
        assert report.ok is False
        assert any('Old BILN' in e for e in report.errors)

    def test_old_biln_validate_ok_with_fmt_biln(self):
        """Old BILN K(1,3) with fmt='biln' → validate returns ok."""
        report = Sequence.validate('K(1,3)-A-K(1,3)', fmt='biln')
        assert report.ok is True


class TestBracketNotation:
    """Tests for .[..] sequential reaction bracket syntax in CABILN.

    Inside .[..], each .Frag(x,y) attaches to the preceding fragment:
      - first entry: host.Rx -> Frag.Ry
      - subsequent entries: prev_frag.Rx -> Frag.Ry
    Auto bond IDs are assigned left-to-right starting from 100.
    """

    # --- _expand_inline_caps unit tests ---

    def test_single_step_bracket_equals_inline_cap(self):
        """Single-step bracket .[boc(4,1)] expands identically to .boc(4,1)."""
        result_bracket, _ = _expand_inline_caps('K.[boc(4,1)]-am')
        result_inline,  _ = _expand_inline_caps('K.boc(4,1)-am')
        assert result_bracket == result_inline

    def test_two_step_bracket_host_annotation(self):
        """Host residue receives only the first bond annotation from a two-step bracket."""
        result, _ = _expand_inline_caps('G.[A(3,1).am(2,1)]')
        # G gets the first auto bond ID at its R3
        assert '(100,3)' in result
        # Second bond ID must NOT appear on G
        assert 'G(100,3)(101' not in result

    def test_two_step_bracket_intermediate_fragment(self):
        """Intermediate fragment receives both incoming and outgoing bond annotations."""
        result, _ = _expand_inline_caps('G.[A(3,1).am(2,1)]')
        # A is intermediate: attached at R1 by bond 100, passes on at R2 by bond 101
        assert 'A(100,1)(101,2)' in result

    def test_two_step_bracket_terminal_fragment(self):
        """Terminal fragment in a two-step bracket receives only its incoming bond."""
        result, _ = _expand_inline_caps('G.[A(3,1).am(2,1)]')
        assert 'am(101,1)' in result

    def test_three_step_bracket_full_chain(self):
        """Three-step bracket produces correct annotations on all three fragments."""
        result, _ = _expand_inline_caps('G.[A(3,1).G(2,1).am(2,1)]')
        assert '(100,3)' in result      # G host → bond 100 at R3
        assert 'A(100,1)(101,2)' in result   # first intermediate
        assert 'G(101,1)(102,2)' in result   # second intermediate
        assert 'am(102,1)' in result         # terminal

    def test_multiple_brackets_get_distinct_bond_ids(self):
        """Two brackets on different residues get non-overlapping auto bond IDs."""
        import re as _re
        result, _ = _expand_inline_caps('C.[trt(4,2)]-G-K.[boc(4,1)]-am')
        ids = [int(i) for i in _re.findall(r'\((\d+),\d+\)', result)
               if int(i) >= 100]
        # Each single-step bracket uses one bond ID; each ID appears exactly twice
        # (once on the host, once on the appended fragment).
        unique_ids = set(ids)
        assert len(unique_ids) == 2, f"Expected 2 distinct auto IDs, got: {unique_ids}"
        # Each ID should appear exactly twice (both bond endpoints)
        for bid in unique_ids:
            assert ids.count(bid) == 2, f"Bond ID {bid} should appear twice: {ids}"

    def test_bracket_and_crosslink_coexist(self):
        """A residue can carry both a .[..] bracket and a .!1 crosslink marker."""
        result, brg = _expand_inline_caps('C.[trt(4,2)].!1(4,4)-G-C.!1')
        # bracket annotation present
        assert '(100,4)' in result
        assert 'trt(100,2)' in result
        # crosslink annotation present
        assert '(!1,4)' in result
        assert brg == {'!1': 4}

    def test_bracket_preserves_other_inline_caps(self):
        """A bracket on one residue does not affect inline caps on other residues."""
        result, _ = _expand_inline_caps('K.[boc(4,1)]-G.ac(2,1)-am')
        # boc bracket: ID 100 on K.R4, boc attached at R1
        assert '(100,4)' in result
        assert 'boc(100,1)' in result
        # ac inline cap: ID 101 on G.R2, ac attached at R1
        assert '(101,2)' in result
        assert 'ac(101,1)' in result

    def test_bracket_in_percent_separated_segment(self):
        """Bracket notation works in a %-separated branch segment."""
        result, _ = _expand_inline_caps('K.!1(4,2)-am%G.[boc(3,1)]-G.!1')
        assert 'boc(' in result
        assert '(!1,4)' in result

    # --- Error cases ---

    def test_empty_bracket_raises(self):
        """Empty .[] raises ValueError."""
        with pytest.raises(ValueError, match='no valid'):
            _expand_inline_caps('C.[]')

    def test_bracket_with_crosslink_marker(self):
        """.[!1(4,4)] is a pure-crosslink bracket — equivalent to inline .!1(4,4)."""
        result, _ = _expand_inline_caps('C.[!1(4,4)]-G-C.!1-am')
        assert '(!1,4)' in result

    def test_bracket_with_monomer_and_crosslinks(self):
        """.[TBMB(4,4).!1(5,4).!2(6,4)] — monomer anchor + crosslink entries."""
        result, brg = _expand_inline_caps(
            'ac-C.[TBMB(4,4).!1(5,4).!2(6,4)]-A-A-C.!1-A-A-C.!2-am')
        assert 'TBMB' in result
        assert '(!1,5)' in result
        assert '(!2,6)' in result
        assert '(!1,4)' in result
        assert '(!2,4)' in result

    def test_bracket_with_trailing_garbage_raises(self):
        """.[trt(4,2)xyz] raises because 'xyz' is not a valid entry."""
        with pytest.raises(ValueError, match='unrecognised'):
            _expand_inline_caps('C.[trt(4,2)xyz]-am')

    def test_bracket_missing_parens_raises(self):
        """.[trt] (no r-group numbers) raises because .trt has no (x,y)."""
        with pytest.raises(ValueError, match='no valid'):
            _expand_inline_caps('C.[trt]-am')

    # --- Backbone peptide branch guard ---

    def test_backbone_r2_r1_pattern_warns_when_two_steps(self):
        """Two or more R2->R1 steps in a bracket emit UserWarning (peptide branch)."""
        with pytest.warns(UserWarning, match='R2->R1|backbone amide|peptide branch'):
            _expand_inline_caps('K.[G(4,1).A(2,1).am(2,1)]')

    def test_single_r2_r1_step_no_warn(self):
        """A single R2->R1 step (e.g. Mal->DBCO style) does NOT warn — ambiguous."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _expand_inline_caps('C.[Mal(4,1).DBCO(2,1)]')

    def test_backbone_pattern_expansion_still_works(self):
        """Warning does not prevent expansion — structure still returned correctly."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            result, _ = _expand_inline_caps('K.[G(4,1).A(2,1).am(2,1)]')
        assert 'G(100,1)(101,2)' in result
        assert 'A(101,1)(102,2)' in result
        assert 'am(102,1)' in result

    def test_non_backbone_r_groups_no_warn(self):
        """R4->R1 (sidechain->cap) inside bracket does NOT warn."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _expand_inline_caps('K.[boc(4,1)]')

    # --- colorize_cabiln ---

    def test_colorize_single_bracket_pair(self):
        """Single bracket pair gets ANSI colour codes on open and close."""
        import re as _re
        coloured = colorize_cabiln('C.[trt(4,2)]')
        # At least one ANSI colour escape is injected
        assert _re.search(r'\033\[\d+m', coloured), "No ANSI escape found"
        # Reset code present after close bracket
        assert '\033[0m' in coloured
        # Original content still intact when escapes stripped
        plain = _re.sub(r'\033\[\d+m', '', coloured)
        assert plain == 'C.[trt(4,2)]'

    def test_colorize_two_pairs_different_colours(self):
        """Two bracket pairs get different ANSI colour codes."""
        coloured = colorize_cabiln('C.[trt(4,2)]-K.[boc(4,1)]')
        # Extract just the escape sequences before each [
        import re as _re
        escapes = _re.findall(r'\033\[\d+m(?=\[)', coloured)
        assert len(escapes) == 2
        assert escapes[0] != escapes[1]

    def test_colorize_use_ansi_false_returns_plain(self):
        """use_ansi=False returns the string unchanged."""
        s = 'C.[trt(4,2)]-K.[boc(4,1)]'
        assert colorize_cabiln(s, use_ansi=False) == s

    def test_colorize_no_brackets_returns_unchanged(self):
        """String with no brackets is returned as-is (modulo no escapes)."""
        s = 'fmoc-C-G-K-am'
        assert colorize_cabiln(s) == s

    # --- End-to-end assembly tests ---

    def test_bracket_single_step_assembly_matches_inline(self):
        """fmoc-C.[trt(4,2)]-am assembles to the same molecule as fmoc-C.trt(4,2)-am."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol_bracket = _romol('fmoc-C.[trt(4,2)]-am')
            mol_inline  = _romol('fmoc-C.trt(4,2)-am')
        from rdkit import Chem
        smi_bracket = Chem.MolToSmiles(mol_bracket)
        smi_inline  = Chem.MolToSmiles(mol_inline)
        assert smi_bracket == smi_inline

    def test_bracket_validate_accepts_syntax(self):
        """Sequence.validate() returns ok=True for a single-step bracket sequence."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            report = Sequence.validate('fmoc-C.[trt(4,2)]-am')
        assert report.ok is True
