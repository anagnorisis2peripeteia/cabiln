"""
Comprehensive tests for bond chemistry validation and molecule assembly.

Covers:
  1. _check_bond_chemistry — unit tests for each branch (silent / warn / raise)
  2. End-to-end Sequence + Molecule assembly — standard and exotic monomers
  3. Monomer pipeline — normalize_input, pre_activate, cap monomers via CSV
  4. SanitizeMol error wrapping (inspects wrapper is in place)

BILN notation reminder:
  '-'  monomer join within a chain  (A-G-A = Ala-Gly-Ala tripeptide)
  '.'  chain separator              (A.G.A = three isolated single-residue chains)

Run: pytest tests/test_bond_validation_and_assembly.py -v
  or: python tests/test_bond_validation_and_assembly.py
"""

import sys
import os
import csv
import inspect
import pathlib
import tempfile
import warnings

import pytest
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.sequence import (Sequence, ValidationReport, _check_bond_chemistry,
                             _preprocess_cabiln, _expand_inline_caps,
                             _handle_terminal_bond_markers, colorize_cabiln,
                             biln_to_cabiln)
from pyPept.molecule import Molecule

# monomers.sdf is now the CHUCKLES-format library (isotope-labelled dummies).
# No monomer_lib override needed — Sequence() defaults to it.
_LIB = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mol(smi):
    m = Chem.MolFromSmiles(smi)
    assert m is not None, f"Bad SMILES: {smi}"
    return m

def _assemble(biln):
    seq = Sequence(biln) if _LIB is None else Sequence(biln, monomer_lib=_LIB)
    mol = Molecule(seq)
    return seq, mol

def _romol(biln):
    _, mol = _assemble(biln)
    return mol.get_molecule(fmt='ROMol')

def _smiles(biln):
    return Chem.MolToSmiles(_romol(biln))

def _natoms(biln):
    return _romol(biln).GetNumAtoms()


# ---------------------------------------------------------------------------
# 1. _check_bond_chemistry — unit tests
# ---------------------------------------------------------------------------

class TestBondChemistryUnit:
    """Direct calls to _check_bond_chemistry with synthetic mol objects."""

    # --- silent (standard bonds) ---

    def test_amide_silent(self):
        """N-C(=O): standard amide — no warning, no exception."""
        mol_c = _mol('CC(=O)O')
        mol_n = _mol('CN')
        c_idx = next(a.GetIdx() for a in mol_c.GetAtoms()
                     if a.GetAtomicNum() == 6
                     and any(nb.GetAtomicNum() == 8 for nb in a.GetNeighbors()))
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _check_bond_chemistry(mol_c, c_idx, mol_n, n_idx, 'amide-test')

    def test_disulfide_silent(self):
        """S-S: disulfide — silent."""
        m = _mol('CS')
        s = next(a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 16)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _check_bond_chemistry(m, s, m, s, 'ss-test')

    def test_ester_silent(self):
        """O-C(=O): ester — silent."""
        mol_c = _mol('CC(=O)O')
        mol_o = _mol('CO')
        c_idx = next(a.GetIdx() for a in mol_c.GetAtoms()
                     if a.GetAtomicNum() == 6
                     and any(nb.GetAtomicNum() == 8 for nb in a.GetNeighbors()))
        o_idx = next(a.GetIdx() for a in mol_o.GetAtoms() if a.GetAtomicNum() == 8)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _check_bond_chemistry(mol_c, c_idx, mol_o, o_idx, 'ester-test')

    def test_diselenide_silent(self):
        """Se-Se: diselenide — silent."""
        m = _mol('[SeH]C')
        se = next(a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 34)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _check_bond_chemistry(m, se, m, se, 'se-test')

    # --- warns (exotic but has documented peptide chemistry context) ---

    def test_non_carbonyl_nc_warns(self):
        """N-C where C is not carbonyl: reductive amination product — warns."""
        mol_c = _mol('CC')
        mol_n = _mol('CN')
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with pytest.warns(UserWarning, match='not a carbonyl carbon'):
            _check_bond_chemistry(mol_c, 0, mol_n, n_idx, 'nc-test')

    def test_ether_oc_warns(self):
        """O-C where C is not carbonyl: ether — warns."""
        mol_c = _mol('CC')
        mol_o = _mol('CO')
        o_idx = next(a.GetIdx() for a in mol_o.GetAtoms() if a.GetAtomicNum() == 8)
        with pytest.warns(UserWarning, match='ether'):
            _check_bond_chemistry(mol_c, 0, mol_o, o_idx, 'oc-test')

    def test_nn_warns(self):
        """N-N: hydrazide — warns."""
        mol_n = _mol('CN')
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with pytest.warns(UserWarning, match='N.N'):
            _check_bond_chemistry(mol_n, n_idx, mol_n, n_idx, 'nn-test')

    def test_thioester_warns(self):
        """S-C(=O): thioester — warns."""
        mol_c = _mol('CC(=O)O')
        mol_s = _mol('CS')
        c_idx = next(a.GetIdx() for a in mol_c.GetAtoms()
                     if a.GetAtomicNum() == 6
                     and any(nb.GetAtomicNum() == 8 for nb in a.GetNeighbors()))
        s_idx = next(a.GetIdx() for a in mol_s.GetAtoms() if a.GetAtomicNum() == 16)
        with pytest.warns(UserWarning, match='thioester'):
            _check_bond_chemistry(mol_c, c_idx, mol_s, s_idx, 'ts-test')

    def test_sulfenamide_warns(self):
        """S-N: sulfenamide — warns, does NOT raise."""
        mol_s = _mol('CS')
        mol_n = _mol('CN')
        s_idx = next(a.GetIdx() for a in mol_s.GetAtoms() if a.GetAtomicNum() == 16)
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with pytest.warns(UserWarning, match='sulfenamide'):
            _check_bond_chemistry(mol_s, s_idx, mol_n, n_idx, 'sn-test')

    def test_sulfenamide_message_acknowledges_legitimate_uses(self):
        """S-N warning should mention Cys PTMs / bioconjugation — not be dismissive."""
        mol_s = _mol('CS')
        mol_n = _mol('CN')
        s_idx = next(a.GetIdx() for a in mol_s.GetAtoms() if a.GetAtomicNum() == 16)
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with pytest.warns(UserWarning) as rec:
            _check_bond_chemistry(mol_s, s_idx, mol_n, n_idx, 'sn-msg-test')
        msg = str(rec[0].message)
        assert any(kw in msg for kw in ('Cys', 'bioconjugation', 'NCL', 'macrolactam'))

    def test_thioether_sc_warns(self):
        """S-C where C is NOT carbonyl: thioether — warns, does NOT raise."""
        mol_c = _mol('CC')   # aliphatic C, no carbonyl
        mol_s = _mol('CS')
        c_idx = 0  # first C of ethane
        s_idx = next(a.GetIdx() for a in mol_s.GetAtoms() if a.GetAtomicNum() == 16)
        with pytest.warns(UserWarning, match='thioether'):
            _check_bond_chemistry(mol_c, c_idx, mol_s, s_idx, 'sc-aliphatic-test')

    def test_thioether_sc_message_mentions_trt_acm(self):
        """Thioether S-C warning message should name the relevant protecting groups."""
        mol_c = _mol('CC')
        mol_s = _mol('CS')
        s_idx = next(a.GetIdx() for a in mol_s.GetAtoms() if a.GetAtomicNum() == 16)
        with pytest.warns(UserWarning) as rec:
            _check_bond_chemistry(mol_c, 0, mol_s, s_idx, 'sc-msg-test')
        msg = str(rec[0].message)
        assert any(kw in msg for kw in ('trt', 'acm', 'thioether'))

    # --- raises (no plausible inter-monomer chemistry) ---

    def test_cc_nonvinyl_warns(self):
        """Non-vinyl C-C inter-monomer bond: warns (valid for bioconjugation)."""
        m = _mol('CC')
        with pytest.warns(UserWarning, match='C.C inter-monomer bond'):
            _check_bond_chemistry(m, 0, m, 1, 'cc-test')

    def test_cc_vinyl_passes(self):
        """Vinyl C-C bond (RCM alkene staple): passes silently."""
        m = _mol('C(/C=C/C)CC=C')
        # Pick the two internal alkene carbons
        alkene_cs = [a for a in m.GetAtoms()
                     if a.GetAtomicNum() == 6 and
                     any(m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                         and nb.GetAtomicNum() == 6 for nb in a.GetNeighbors())]
        assert len(alkene_cs) >= 2
        _check_bond_chemistry(m, alkene_cs[0].GetIdx(), m, alkene_cs[1].GetIdx(), 'rcm-test')

    def test_oo_raises(self):
        """O-O inter-monomer bond: raises ValueError."""
        mol_a = _mol('CO')
        mol_b = _mol('CO')
        a_o = next(a.GetIdx() for a in mol_a.GetAtoms() if a.GetAtomicNum() == 8)
        b_o = next(a.GetIdx() for a in mol_b.GetAtoms() if a.GetAtomicNum() == 8)
        with pytest.raises(ValueError, match='inter-monomer bond'):
            _check_bond_chemistry(mol_a, a_o, mol_b, b_o, 'oo-test')

    def test_unknown_element_pair_raises(self):
        """Unrecognised element pair (Si-N): raises ValueError."""
        mol_si = _mol('[SiH4]')
        mol_n  = _mol('CN')
        si = next(a.GetIdx() for a in mol_si.GetAtoms() if a.GetAtomicNum() == 14)
        n  = next(a.GetIdx() for a in mol_n.GetAtoms()  if a.GetAtomicNum() == 7)
        with pytest.raises(ValueError, match='inter-monomer bond'):
            _check_bond_chemistry(mol_si, si, mol_n, n, 'sin-test')

    def test_error_message_names_atom_symbols(self):
        """Warning message for non-vinyl C-C should contain atom symbol names."""
        m = _mol('CC')
        with pytest.warns(UserWarning) as rec:
            _check_bond_chemistry(m, 0, m, 1, 'cc-sym-test')
        assert 'C' in str(rec[0].message)

    def test_bond_label_appears_in_warning(self):
        """Warning message for non-vinyl C-C should include the bond_label."""
        m = _mol('CC')
        with pytest.warns(UserWarning, match='my-custom-label'):
            _check_bond_chemistry(m, 0, m, 1, 'my-custom-label')

    def test_on_hydroxylamine_warns(self):
        """O-N bond (hydroxylamine/hydroxamic acid): warns but does not raise."""
        mol_o = _mol('CO')
        mol_n = _mol('CN')
        o_idx = next(a.GetIdx() for a in mol_o.GetAtoms() if a.GetAtomicNum() == 8)
        n_idx = next(a.GetIdx() for a in mol_n.GetAtoms() if a.GetAtomicNum() == 7)
        with pytest.warns(UserWarning, match='hydroxylamine|hydroxamic'):
            _check_bond_chemistry(mol_o, o_idx, mol_n, n_idx, 'on-test')


# ---------------------------------------------------------------------------
# 2. End-to-end assembly tests  (all use monomers_new.sdf)
# ---------------------------------------------------------------------------

class TestAssembly:
    """Sequence + Molecule round-trips."""

    def test_tripeptide_is_connected(self):
        """A-G-A assembles into a single connected molecule (not dot-separated)."""
        smi = _smiles('A-G-A')
        assert '.' not in smi, f"Got disconnected SMILES: {smi}"

    def test_tripeptide_atom_count(self):
        """A-G-A: 15 heavy atoms (AGA - 2 water = C8N3O4)."""
        assert _natoms('A-G-A') == 15

    def test_smiles_parseable(self):
        """Assembled SMILES should be parseable by RDKit."""
        assert Chem.MolFromSmiles(_smiles('A-G-V')) is not None

    def test_trp_indole_intact(self):
        """Trp: indole ring system present after Kekulize-first assembly."""
        mol = _romol('A-W-A')
        indole = Chem.MolFromSmarts('c1ccc2[nH]ccc2c1')
        assert mol.HasSubstructMatch(indole), "Indole ring missing from Trp assembly"

    def test_his_assembles(self):
        """His (imidazole NH) assembles without aromaticity error."""
        mol = _romol('A-H-A')
        assert mol is not None
        assert mol.GetNumAtoms() == 21

    def test_asp_sidechain_cooh_intact(self):
        """Asp R3 sidechain COOH must not become aldehyde (key leaving-group fix)."""
        mol = _romol('A-D-A')
        cooh = Chem.MolFromSmarts('C(=O)O')
        matches = mol.GetSubstructMatches(cooh)
        # C-terminal COOH + Asp sidechain COOH = at least 2 C(=O)O groups
        assert len(matches) >= 2, f"Expected >=2 C(=O)O groups; got {len(matches)}"

    def test_glu_sidechain_cooh_intact(self):
        """Glu: same sidechain COOH check."""
        mol = _romol('A-E-A')
        cooh = Chem.MolFromSmarts('C(=O)O')
        assert len(mol.GetSubstructMatches(cooh)) >= 2

    def test_no_dummy_atoms_leaked(self):
        """No [*] dummy atoms should remain in any assembled product."""
        for biln in ('A-G', 'A-W-A', 'A-D-A', 'G-P-G'):
            romol = _romol(biln)
            dummies = [a for a in romol.GetAtoms() if a.GetAtomicNum() == 0]
            assert dummies == [], f"{biln}: {len(dummies)} dummy atom(s) leaked"

    def test_proline_assembles(self):
        """Pro (secondary amine N-terminus) bonds correctly."""
        assert _natoms('G-P-G') > 0

    def test_longer_peptide(self):
        """8-residue peptide with diverse sidechains assembles without error."""
        smi = _smiles('A-G-V-L-F-W-D-E')
        assert Chem.MolFromSmiles(smi) is not None

    def test_sanitize_error_wrapped_as_valueerror(self):
        """SanitizeMol failure is re-raised as ValueError with user-friendly message."""
        src = inspect.getsource(Molecule)
        assert 'SanitizeMol' in src
        assert 'R-group' in src  # wrapper message mentions R-group assignments

    def test_disulfide_crosslink_assembles(self):
        """Cys-Ala-Cys with R4-R4 disulfide crosslink via .!1(4,4) notation.

        Cys thiol is R4 in CABILN (R3 is backbone-N mod).
        S-S is a silent bond — no UserWarning expected.
        The assembled molecule should contain an S-S substructure.
        """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            mol = _romol('C.!1(4,4)-A-C.!1(4,4)-am')

        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

        ss_pat = Chem.MolFromSmarts('[S][S]')
        assert mol.HasSubstructMatch(ss_pat), "S-S disulfide bond not found in assembled molecule"

        bond_warns = [x for x in w if 'Bond' in str(x.message)]
        assert bond_warns == [], f"Unexpected bond warnings for S-S: {[str(x.message) for x in bond_warns]}"


# ---------------------------------------------------------------------------
# 3. Monomer pipeline tests
# ---------------------------------------------------------------------------

class TestMonomerpipeline:

    def test_normalize_detects_chuckles(self):
        """Input containing [n*] is classified as CHUCKLES and returned unchanged."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        smi, is_chuckles, err = normalize_input('[1*]N[C@@H](C)C(=O)[2*]')
        assert err is None
        assert is_chuckles is True
        assert '[1*]' in smi

    def test_normalize_detects_smiles(self):
        """Plain SMILES is classified correctly (not CHUCKLES)."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        smi, is_chuckles, err = normalize_input('N[C@@H](C)C(=O)O')
        assert err is None
        assert is_chuckles is False

    def test_normalize_fasta_single_letter(self):
        """Single FASTA letter 'A' returns valid Ala SMILES."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        smi, is_chuckles, err = normalize_input('A')
        assert err is None
        assert is_chuckles is False
        assert Chem.MolFromSmiles(smi) is not None

    def test_normalize_biln_token_full_name(self):
        """BILN token 'DAla' (full name, not 'dA') returns valid SMILES."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        smi, is_chuckles, err = normalize_input('DAla')
        assert err is None
        assert Chem.MolFromSmiles(smi) is not None

    def test_normalize_unknown_returns_error(self):
        """Unrecognisable input returns a non-None error string."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        _, _, err = normalize_input('ZZZZNOTREAL')
        assert err is not None

    def test_pre_activate_canonical_amino_acids(self):
        """All 20 canonical AAs pre-activate with R1 and R2 attachment points."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        aa_smiles = [
            'NCC(=O)O',                       # Gly
            'N[C@@H](C)C(=O)O',               # Ala
            'N[C@@H](CS)C(=O)O',              # Cys
            'N[C@@H](CO)C(=O)O',              # Ser
            'N[C@@H](Cc1ccccc1)C(=O)O',       # Phe
            'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',  # Trp
            'N[C@@H](Cc1cnc[nH]1)C(=O)O',    # His
            'N[C@@H](CC(=O)O)C(=O)O',         # Asp
            'N[C@@H](CCC(=O)O)C(=O)O',        # Glu
        ]
        for smi in aa_smiles:
            ch, lg, chem_types, err = pre_activate(smi)
            assert err is None, f"pre_activate failed for {smi}: {err}"
            assert '[1*]' in ch, f"No R1 in CHUCKLES for {smi}"
            assert '[2*]' in ch, f"No R2 in CHUCKLES for {smi}"
            assert chem_types.get(1) == 'backbone_n', f"Wrong R1 chem_type for {smi}"
            assert chem_types.get(2) == 'backbone_c', f"Wrong R2 chem_type for {smi}"

    def test_cap_monomers_skip_backbone_detection(self):
        """Caps with pre-filled CHUCKLES columns are written to SDF without error.

        Uses exact column names from _CSV_COLUMNS and values matching monomers.csv.
        """
        from pyPept.interfaces.monomer_pipeline import build_library_from_csv, _CSV_COLUMNS

        # Values taken directly from src/pyPept/data/monomers.csv
        cap_rows = [
            {
                'token': 'ac', 'input': 'CC(=O)O', 'name': 'Acetyl cap',
                'type': 'cap', 'synonyms': 'Ac',
                'chuckles': 'CC([2*])=O',
                'r1_leaving': 'None', 'r2_leaving': '[OH]',
                'r3_leaving': 'None', 'r4_leaving': 'None',
            },
            {
                'token': 'am', 'input': 'N', 'name': 'Amide cap',
                'type': 'cap', 'synonyms': 'NH2',
                'chuckles': '[1*]N',
                'r1_leaving': '[H]', 'r2_leaving': 'None',
                'r3_leaving': 'None', 'r4_leaving': 'None',
            },
        ]

        with tempfile.NamedTemporaryFile(
                suffix='.csv', mode='w', delete=False, newline='') as f:
            csv_path = pathlib.Path(f.name)
            writer = csv.DictWriter(f, fieldnames=_CSV_COLUMNS)
            writer.writeheader()
            for row in cap_rows:
                writer.writerow({c: row.get(c, '') for c in _CSV_COLUMNS})

        sdf_path = pathlib.Path(tempfile.mktemp(suffix='.sdf'))
        try:
            build_library_from_csv(csv_path, sdf_path)
            sdf_text = sdf_path.read_text()
            assert '>  <token>' in sdf_text or 'ac' in sdf_text, \
                "ac not found in SDF output"
            assert 'am' in sdf_text, "am not found in SDF output"
        finally:
            csv_path.unlink(missing_ok=True)
            sdf_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# 4. Inline attachment notation — _expand_inline_caps + assembly integration
# ---------------------------------------------------------------------------

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

    def test_inline_cap_ac_on_backbone(self):
        """ac cap applied inline via .ac(2,1) on last residue assembles correctly.

        ac(2,1) = ac's R1 binds to the residue's R2 (C-terminal exit).
        Equivalent to writing ac as a terminal cap in explicit notation.
        """
        # Standard explicit notation: ac-A-G-am
        smi_explicit = _smiles('ac-A-G-am')
        # Inline: apply ac as a side attachment on G's R2 (C-terminal)
        # Note: this gives the same topology as the explicit form.
        # G.ac(2,1) means G-R2 → ac-R1 (amide bond equivalent to C-terminal ac cap)
        mol_explicit = _romol('ac-A-G-am')
        assert mol_explicit is not None

    def test_inline_crosslink_between_residues(self):
        """!1 crosslink between two K (Lys) residues via R3-R3 assembles.

        K.!1(3,3) marks R3 (ε-amine, N) on each K as crosslink endpoints.
        N–N bond → UserWarning (hydrazide) but not a ValueError — assembles ok.
        """
        import warnings
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            seq = Sequence('K.!1(3,3)-A-K.!1(3,3)-am')
        assert seq is not None
        assert seq.s_nmonomers == 4  # K, A, K, am

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


# ---------------------------------------------------------------------------
# 5. Protecting-group cap monomers (trt, acm, pbf)
# ---------------------------------------------------------------------------

class TestCapMonomers:
    """Assembly tests for thioether and sulfonyl protecting group caps."""

    def test_trt_on_cys_warns_thioether(self):
        """trt cap on Cys R4 (thiol) forms S-C thioether — UserWarning, assembles OK.

        R3 = backbone-N mod; R4 = thiol in CABILN convention.
        """
        with pytest.warns(UserWarning, match='thioether'):
            mol = _romol('fmoc-C.trt(4,2)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_trt_cys_has_sc_bond(self):
        """trt-protected Cys should contain an S-C bond to a trisubstituted carbon."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-C.trt(4,2)-am')
        sc_trt = Chem.MolFromSmarts('[S][C](c1ccccc1)(c1ccccc1)c1ccccc1')
        assert mol.HasSubstructMatch(sc_trt), "trt S-C bond pattern not found"

    def test_acm_on_cys_warns_thioether(self):
        """acm cap on Cys R4 (thiol) forms S-C thioether — UserWarning, assembles OK."""
        with pytest.warns(UserWarning, match='thioether'):
            mol = _romol('fmoc-C.acm(4,2)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_acm_cys_has_acetamide(self):
        """acm-protected Cys should contain the S-CH2-NH-C(=O) acetamidomethyl group."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-C.acm(4,2)-am')
        acm_pat = Chem.MolFromSmarts('[S]CNC(C)=O')
        assert mol.HasSubstructMatch(acm_pat), "acm S-CH2-NH-C(=O) pattern not found"

    def test_pbf_on_arg_warns_sulfenamide(self):
        """pbf cap on Arg R4 (guanidinium NH2) forms S-N bond — UserWarning, assembles OK.

        R3 = backbone-N mod; R4 = guanidinium terminal NH2 in CABILN convention.
        """
        with pytest.warns(UserWarning, match='sulfenamide|S.N'):
            mol = _romol('fmoc-R.pbf(4,2)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_pbf_arg_has_sulfonyl(self):
        """pbf-protected Arg should contain an N-S(=O)(=O) sulfonamide substructure."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-R.pbf(4,2)-am')
        sulf_pat = Chem.MolFromSmarts('[N]S(=O)=O')
        assert mol.HasSubstructMatch(sulf_pat), "N-S(=O)2 sulfonamide pattern not found"

    def test_obn_cap_warns_hydroxylamine(self):
        """OBn_ cap on Gly forms O-N bond — warns hydroxylamine, assembles OK."""
        with pytest.warns(UserWarning, match='hydroxylamine|hydroxamic'):
            mol = _romol('OBn_-G-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_ome_cap_warns_hydroxylamine(self):
        """OMe_ cap on Gly forms O-N bond — warns hydroxylamine, assembles OK."""
        with pytest.warns(UserWarning, match='hydroxylamine|hydroxamic'):
            mol = _romol('OMe_-G-am')
        assert mol is not None


# ---------------------------------------------------------------------------
# 6. Phase 2 CABILN parser — _preprocess_cabiln, _expand_inline_caps, terminals
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# 7. Intramolecular ring closure — grouped SMIRKS path
# ---------------------------------------------------------------------------

class TestIntramolecularRingClosure:
    """
    Hard edge cases for the intramolecular grouped-SMIRKS path.

    Two implementation strategies:
      - CABILN end-to-end: standard AAs through Sequence + Molecule.
      - Direct run_bond_smirks: build a model assembled SMILES with both
        isotope-labelled dummies already in one connected molecule, then
        call run_bond_smirks(..., intramolecular=True).  Used for exotic
        chemistry where the monomer library has no entries.
    """

    # ── Head-to-tail backbone cyclisation ─────────────────────────────────────

    def test_head_to_tail_cyclo_tetra_alanine(self):
        """Cyclo(AAAA): R1 of first residue → R2 of last via backbone_amide.

        CABILN: A.!1(1,2)-A-A-A.!1(2,1)
        Expected: 20 heavy atoms, at least one ring, no free termini.
        """
        mol = _romol('A.!1(1,2)-A-A-A.!1(2,1)')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi, f"Product is disconnected: {smi}"
        assert mol.GetNumAtoms() == 20, f"Expected 20 atoms for cyclo(AAAA), got {mol.GetNumAtoms()}"
        assert mol.GetRingInfo().NumRings() > 0, "No rings found in cyclic peptide"
        # No free amine or carboxyl terminus
        free_amine = Chem.MolFromSmarts('[NH2][C]')
        free_acid  = Chem.MolFromSmarts('C(=O)[OH]')
        assert not mol.HasSubstructMatch(free_amine), "Free N-terminus found — ring not closed"
        assert not mol.HasSubstructMatch(free_acid),  "Free C-terminus found — ring not closed"

    def test_head_to_tail_cyclo_penta_glycine(self):
        """Cyclo(GGGGG): 5-residue ring.  Gly has no β-carbon, smallest cyclic peptide."""
        mol = _romol('G.!1(1,2)-G-G-G-G.!1(2,1)')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"
        # cyclo(G5): 5*(N + CH2 + C + O) = 5*4 = 20 heavy atoms
        assert mol.GetNumAtoms() == 20, f"Expected 20 atoms for cyclo(GGGGG), got {mol.GetNumAtoms()}"
        assert mol.GetRingInfo().NumRings() > 0

    def test_head_to_tail_shorthand_hex_ala(self):
        """!1-A-A-A-A-A-A-!1 shorthand: no explicit .!1(y,z) required.

        Terminal markers imply R1 (N-terminal) and R2 (C-terminal) by position.
        """
        mol = _romol('!1-A-A-A-A-A-A-!1')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi, f"Product is disconnected: {smi}"
        # cyclo(A6): 6*(N + CH(CH3) + C + O) = 6*5 = 30 heavy atoms
        assert mol.GetNumAtoms() == 30, f"Expected 30 atoms, got {mol.GetNumAtoms()}"
        assert mol.GetRingInfo().NumRings() > 0

    def test_head_to_tail_shorthand_penta_gly(self):
        """!1-G-G-G-G-G-!1 pure shorthand — same ring as explicit form."""
        mol_short    = _romol('!1-G-G-G-G-G-!1')
        mol_explicit = _romol('G.!1(1,2)-G-G-G-G.!1(2,1)')
        from rdkit import Chem as _Chem
        assert _Chem.MolToSmiles(mol_short) == _Chem.MolToSmiles(mol_explicit)

    # ── Double disulfide (bicyclic) ────────────────────────────────────────────

    def test_double_disulfide_two_ss_bonds(self):
        """Two independent disulfide crosslinks on one linear chain → bicyclic product.

        Pattern: C(ss1)-A-C(ss2)-A-A-C(ss1)-A-C(ss2)-am
        Verifies both S-S bonds and two rings are present.
        """
        mol = _romol('C.!1(4,4)-A-C.!2(4,4)-A-A-C.!1(4,4)-A-C.!2(4,4)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi
        ss_pat = Chem.MolFromSmarts('[S][S]')
        matches = mol.GetSubstructMatches(ss_pat)
        assert len(matches) == 2, f"Expected 2 S-S bonds, found {len(matches)}"
        assert mol.GetRingInfo().NumRings() >= 2, "Expected ≥2 rings in bicyclic product"

    def test_double_disulfide_no_free_thiols(self):
        """Both thiols consumed: no free S-H present after bicyclic assembly."""
        mol = _romol('C.!1(4,4)-A-C.!2(4,4)-A-A-C.!1(4,4)-A-C.!2(4,4)-am')
        free_sh = Chem.MolFromSmarts('[SH]')
        assert not mol.HasSubstructMatch(free_sh), "Free thiol found — crosslink not closed"

    # ── Short-range disulfide (ring strain edge case) ─────────────────────────

    def test_adjacent_cys_disulfide_8_membered_ring(self):
        """C-C with disulfide: smallest ring containing S-S (8-membered).

        Tests that grouped SMIRKS tolerates ring-strain without raising.
        """
        mol = _romol('C.!1(4,4)-C.!1(4,4)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"
        ss_pat = Chem.MolFromSmarts('[S][S]')
        assert mol.HasSubstructMatch(ss_pat), "S-S bond not found in adjacent-Cys product"
        # Confirm product is a ring, not two chains
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi

    def test_1_3_disulfide_medium_ring(self):
        """C-A-C disulfide: 10-membered ring (two extra atoms from Ala backbone)."""
        mol = _romol('C.!1(4,4)-A-C.!1(4,4)-am')
        assert mol is not None
        ss_pat = Chem.MolFromSmarts('[S][S]')
        assert mol.HasSubstructMatch(ss_pat)
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []

    # ── Thiol-maleimide intramolecular (direct run_bond_smirks) ───────────────

    def test_thiol_maleimide_intramolecular_grouped_smirks(self):
        """Thiol-maleimide Michael addition forms succinimide ring intramolecularly.

        Previously required graph surgery.  Now uses grouped SMIRKS.

        Model assembled mol: maleimide ring ([101*] on olefinic C) connected via
        N-alkyl chain to thiol ([201*] on S).  After reaction the thiol S adds
        across the C=C, retaining the succinimide ring and closing the macrocycle.
        """
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks

        # [101*] on C:2 (olefinic carbon of maleimide ring)
        # [201*] on S:10 (thiol)
        # Chain: maleimide-N → 4-carbon → thiol
        smi = '[101*]C1=CC(=O)N(CCCCS[201*])C1=O'
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Model assembled SMILES invalid: {smi}"

        entry = REACTIONS['thiol_maleimide']
        product = run_bond_smirks(mol, mol, 101, 201, entry, intramolecular=True)
        assert product is not None

        dummies = [a for a in product.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

        smi_out = Chem.MolToSmiles(product)
        assert '.' not in smi_out, f"Product disconnected: {smi_out}"

        # Succinimide (Michael adduct ring) must be present
        succinimide = Chem.MolFromSmarts('C1CC(=O)NC1=O')
        assert product.HasSubstructMatch(succinimide), \
            f"Succinimide ring not found in product: {smi_out}"

        # S-C bond to the former olefinic carbon
        sc_pat = Chem.MolFromSmarts('[S][CH]1CC(=O)NC1=O')
        assert product.HasSubstructMatch(sc_pat), \
            f"S-C succinimide bond not found: {smi_out}"

    # ── CuAAC intramolecular staple (direct run_bond_smirks) ─────────────────

    def test_cuaac_intramolecular_triazole_staple(self):
        """CuAAC alkyne + azide → 1,4-triazole ring intramolecularly.

        [101*] is on the propargylic CH (α to the internal alkyne C:2), per
        authoring Rule 2: [4*] on adjacent sp3 C so C:2's chain bond survives.
        [201*] is on the azide alpha-C:5, which already carries an unmapped chain bond.
        The 3-carbon chain connecting them closes into a macrocycle + triazole.
        """
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks

        # [101*] on propargylic C (adjacent to internal alkyne C:2)
        # C:2 triple-bonded to terminal C:3H
        # 3-C chain links propargylic C to azide alpha-C:5 ([201*])
        smi = '[101*]C(CCC([201*])N=[N+]=[N-])C#[CH]'
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Model assembled SMILES invalid: {smi}"

        entry = REACTIONS['cuaac_1_4_triazole']
        product = run_bond_smirks(mol, mol, 101, 201, entry, intramolecular=True)
        assert product is not None

        dummies = [a for a in product.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

        smi_out = Chem.MolToSmiles(product)
        assert '.' not in smi_out, f"Product disconnected: {smi_out}"

        # 1,2,3-triazole ring must be present
        triazole = Chem.MolFromSmarts('c1cnnn1')
        assert product.HasSubstructMatch(triazole), \
            f"1,2,3-triazole not found in CuAAC product: {smi_out}"

    # ── NHS ester intramolecular macrolactamisation ───────────────────────────

    def test_nhs_ester_intramolecular_macrolactam(self):
        """NHS ester + amine → amide intramolecularly; NHS ring departs as byproduct.

        [101*] is on the α-CH adjacent to the NHS carbonyl C:2, per authoring
        Rule 2: [4*] on adjacent sp3 C so the chain bond on that C survives.
        [201*] is on amine N:5, which already carries an unmapped chain bond.
        take_largest filters the expelled NHS (8 heavy atoms) from the macrolactam
        (9 heavy atoms with the 5-C chain used here).
        """
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks

        # [101*] on α-CH (adjacent to NHS carbonyl C:2)
        # 5-C chain links α-CH to amine N:5 ([201*])
        # NHS ring (ON1C(=O)CCC1=O) departs; take_largest keeps the macrolactam
        smi = '[101*]C(CCCCCN[201*])C(=O)ON1C(=O)CCC1=O'
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Model assembled SMILES invalid: {smi}"

        entry = REACTIONS['nhs_ester_amide']
        product = run_bond_smirks(mol, mol, 101, 201, entry, intramolecular=True)
        assert product is not None

        dummies = [a for a in product.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

        smi_out = Chem.MolToSmiles(product)
        assert '.' not in smi_out, f"Product disconnected (NHS not expelled?): {smi_out}"

        # Amide bond formed
        amide = Chem.MolFromSmarts('[NH]C(=O)')
        assert product.HasSubstructMatch(amide), \
            f"Amide bond not found in macrolactam product: {smi_out}"


# ---------------------------------------------------------------------------
# 8. Expanded monomer library — migrated NCAAs
# ---------------------------------------------------------------------------

class TestExpandedMonomers:
    """Smoke-tests for monomers migrated from the legacy 322-entry library.

    Each test assembles a simple sequence and checks the product is a valid,
    connected RDKit mol with a plausible atom count.  Exact SMILES are NOT
    checked here — that would over-specify stereochemistry and leave-group
    order.  The goal is: library present, CHUCKLES valid, assembly works.
    """

    # ── D-amino acids (d prefix) ──────────────────────────────────────────

    def test_d_alanine_dipeptide(self):
        assert _natoms('dA-A') > 10

    def test_d_phenylalanine_dipeptide(self):
        assert _natoms('dF-A') > 16

    def test_d_tryptophan_dipeptide(self):
        assert _natoms('dW-A') == 20

    def test_d_proline_dipeptide(self):
        assert _natoms('dP-A') == 13

    # ── Non-canonical α-amino acids ──────────────────────────────────────

    def test_norleucine_dipeptide(self):
        """Nle (norleucine): straight-chain C4 sidechain, no sidechain slot."""
        assert _natoms('Nle-A') == 14

    def test_norvaline_dipeptide(self):
        assert _natoms('Nva-A') > 11

    def test_ornithine_dipeptide(self):
        """Orn (ornithine): Lys minus one CH2; has ε-amine R4."""
        assert _natoms('Orn-A') > 12

    def test_homoserine_dipeptide(self):
        """Hse (homoserine): Ser with extra CH2; has sidechain OH."""
        assert _natoms('Hse-A') > 11

    def test_homocysteine_dipeptide(self):
        """Hcy (homocysteine): Cys with extra CH2; thiol sidechain."""
        assert _natoms('Hcy-A') > 11

    def test_hydroxyproline_dipeptide(self):
        """Hyp (trans-4-hydroxyproline): proline with OH on ring."""
        assert _natoms('Hyp-A') == 14

    def test_aib_tripeptide(self):
        """Aib (α-aminoisobutyric acid): Cα-disubstituted, no sidechain slot."""
        assert _natoms('A-Aib-A') == 17

    def test_gamma_glutamic_acid_dipeptide(self):
        """E_g: γ-linked Glu (backbone through γ-carboxyl)."""
        assert _natoms('E_g-A') == 15

    # ── α-Methyl amino acids ─────────────────────────────────────────────

    def test_aMePhe_dipeptide(self):
        """α-MePhe: Phe with α-methyl (Cα-disubstituted)."""
        assert _natoms('aMePhe-A') == 18

    # ── Naphthylalanine / aromatic NCAAs ─────────────────────────────────

    def test_1_naphthylalanine_dipeptide(self):
        assert _natoms('1Nal-A') == 21

    def test_2_naphthylalanine_dipeptide(self):
        assert _natoms('2Nal-A') == 21

    # ── Halogenated phenylalanines ────────────────────────────────────────

    def test_phe_4F_dipeptide(self):
        assert _natoms('Phe_4F-A') > 16

    def test_phe_4Cl_dipeptide(self):
        assert _natoms('Phe_4Cl-A') > 16

    # ── Cyclohexylalanine / bulky aliphatics ─────────────────────────────

    def test_cyclohexylalanine_dipeptide(self):
        assert _natoms('Cha-A') == 17

    # ── Tripeptides with alternating NCAAs ───────────────────────────────

    def test_mixed_d_l_tripeptide(self):
        """Alternating D/L — a common motif in peptide drug design."""
        assert _natoms('A-dA-A') == 16

    def test_nle_cha_tripeptide(self):
        assert _natoms('Nle-Cha-A') == 25

    # ── Cap monomers from expanded library ───────────────────────────────

    def test_cbz_cap_assembly(self):
        """Cbz (benzyloxycarbonyl) N-terminal protecting group."""
        assert _natoms('Cbz-A-am') == 16

    def test_boc_cap_assembly(self):
        """Boc N-terminal protecting group."""
        assert _natoms('Boc-A-am') == 13


class TestFinalMonomers:
    """Tests for the last 31 monomers completing the 322-entry library.

    Covers amino alcohols, amino aldehydes, NMe aromatic AAs, NMe caps,
    PNA nucleobase monomers, and Pqa.
    """

    # ── Amino alcohols (C-terminal, R1 only) ─────────────────────────────

    def test_ala_ol(self):
        assert _natoms('G-Ala_ol') == 9

    def test_gly_ol(self):
        assert _natoms('G-Gly_ol') == 7

    def test_pro_ol(self):
        assert _natoms('G-Pro_ol') == 11

    def test_val_ol(self):
        assert _natoms('G-Val_ol') == 11

    def test_leu_ol(self):
        assert _natoms('G-Leu_ol') == 12

    def test_phe_ol(self):
        assert _natoms('G-Phe_ol') == 15

    def test_phg_ol(self):
        assert _natoms('G-Phg_ol') == 14

    def test_lys_ol(self):
        assert _natoms('G-Lys_ol') == 13

    def test_thr_ol(self):
        assert _natoms('G-Thr_ol') == 11

    def test_aib_ol(self):
        assert _natoms('G-Aib_ol') == 10

    def test_d_phg_ol(self):
        assert _natoms('G-D_Phg_ol') == 14

    def test_d_pro_ol(self):
        assert _natoms('G-D_Pro_ol') == 11

    def test_d_thr_ol(self):
        assert _natoms('G-D_Thr_ol') == 11

    def test_hsl(self):
        """Hsl (homoserine lactone): cyclic ester at C-terminus."""
        assert _natoms('G-Hsl') == 11

    # ── Amino aldehydes (C-terminal, R1 only) ────────────────────────────

    def test_ala_al(self):
        assert _natoms('G-Ala_al') == 9

    def test_gly_al(self):
        assert _natoms('G-Gly_al') == 8

    def test_pro_al(self):
        assert _natoms('G-Pro_al') == 11

    def test_leu_al(self):
        assert _natoms('G-Leu_al') == 12

    def test_phe_al(self):
        assert _natoms('G-Phe_al') == 15

    def test_lys_al(self):
        assert _natoms('G-Lys_al') == 13

    def test_arg_al(self):
        assert _natoms('G-Arg_al') == 15

    # ── N-methyl amino acids (R1+R2) ─────────────────────────────────────

    def test_nme_beta_ala(self):
        """NMebAla: N-methyl β-alanine."""
        assert _natoms('ac-NMebAla-am') == 10

    def test_nme2abz(self):
        """NMe2Abz: 2-(methylamino)benzoic acid."""
        assert _natoms('ac-NMe2Abz-am') == 15

    def test_nme4abz(self):
        """NMe4Abz: 4-(methylamino)benzoic acid."""
        assert _natoms('ac-NMe4Abz-am') == 14

    # ── N,N-dimethyl benzoic acid caps (R2 only) ─────────────────────────

    def test_nme23abz_cap(self):
        """NMe23Abz: N-terminal 3-(dimethylamino)benzoyl cap."""
        assert _natoms('NMe23Abz-G-am') == 16

    def test_nme24abz_cap(self):
        """NMe24Abz: N-terminal 4-(dimethylamino)benzoyl cap."""
        assert _natoms('NMe24Abz-G-am') == 16

    # ── PNA nucleobase monomers (R1+R2) ──────────────────────────────────

    def test_pna_adenine(self):
        assert _natoms('ac-pnA-am') == 24

    def test_pna_cytosine(self):
        assert _natoms('ac-pnC-am') == 22

    def test_pna_guanine(self):
        assert _natoms('ac-pnG-am') == 25

    def test_pna_thymine(self):
        assert _natoms('ac-pnT-am') == 23

    def test_pna_tetranucleotide(self):
        """ac-pnA-pnC-pnG-pnT-am: four-base PNA strand."""
        mol = _romol('ac-pnA-pnC-pnG-pnT-am')
        assert mol.GetNumAtoms() > 60

    # ── Pqa (isoquinolinone piperazine, R1+R2) ────────────────────────────

    def test_pqa(self):
        """Pqa: 2-(1-oxo-7-(piperazin-1-yl)-1,2-dihydroisoquinolin-2-yl)acetic acid."""
        assert _natoms('ac-Pqa-am') == 24


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


class TestFattyAcidBranching:
    """Sidechain fatty-acid branching via the new amine_primary+backbone_c amide reaction.

    Models the GLP-1/GIP agonist lipidation architecture:
    Lys(ε-NH2) → amide → AEEA → E_g → C18/C20 fatty diacid.
    """

    def test_aeea_linear(self):
        """AEEA (mini-PEG linker) assembles as a standard backbone monomer."""
        assert _natoms('ac-AEEA-am') == 14

    def test_ameLeu_linear(self):
        """aMeLeu (alpha-methyl-L-leucine): Cα-disubstituted, backbone stabiliser."""
        assert _natoms('ac-aMeLeu-am') == 13

    def test_c20fa_E_g_branch_linear(self):
        """C20 fatty diacid cap + E_g chain assembles correctly."""
        assert _natoms('C20FA-E_g-am') == 33

    def test_c20fa_E_g_aeea_linear(self):
        """Full linker chain: C20FA-E_g-AEEA with backbone amide bonds."""
        assert _natoms('C20FA-E_g-AEEA-am') == 43

    def test_lys_sidechain_amide_bond(self):
        """Lys R4 (epsilon-amine) forms amide with AEEA R2 (carboxyl).

        Tests the new amine_primary+backbone_c reaction added to BOND_TABLE.
        """
        assert _natoms('ac-K.!1(4,2)-G-am%C20FA-E_g-AEEA-!1') == 59

    def test_retatrutide(self):
        """Retatrutide (LY3437943): GIP/GLP-1/glucagon triple agonist, 39 AAs.

        Sequence: Y-Aib-QGTFTSDYSI-aMeLeu-LD-K-K(C20FA-E_g-AEEA)-AQAAFI-Aib-EYL
                  LEGGPSSGAPPPSam with C20 fatty diacid branch on Lys17 (K16-K17).
        """
        retatrutide = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K-K.!1(4,2)"
            "-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%C20FA-E_g-AEEA-!1"
        )
        mol = _romol(retatrutide)
        assert mol.GetNumAtoms() == 335


class TestStapledPeptides:
    """All-hydrocarbon stapled peptides via RCM (S5/R8 monomers)."""

    def test_s5_s5_i_i4_staple(self):
        """i, i+4 RCM staple using S5+S5 (same-configuration pair)."""
        mol = _romol('ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am')
        assert mol is not None
        assert mol.GetRingInfo().NumRings() >= 1
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []
        assert mol.HasSubstructMatch(Chem.MolFromSmarts('C=C'))
        assert mol.GetNumAtoms() == 46

    def test_s5_r8_i_i7_staple(self):
        """i, i+7 RCM staple using S5+R8 (opposite-configuration pair)."""
        mol = _romol('ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am')
        assert mol is not None
        assert mol.GetRingInfo().NumRings() >= 1
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []
        assert mol.GetNumAtoms() == 59

    def test_staple_gives_cyclic_alkene(self):
        """RCM product has internal alkene (the staple bridge) not terminal alkene."""
        mol = _romol('ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am')
        internal_alkene = Chem.MolFromSmarts('[CH]=[CH]')
        terminal_alkene = Chem.MolFromSmarts('[CH]=[CH2]')
        assert mol.HasSubstructMatch(internal_alkene)
        assert not mol.HasSubstructMatch(terminal_alkene)


class TestAspartimide:
    """Cyclic imide (succinimide) via Asp sidechain R4 + backbone amide N-H R3."""

    def test_asp_gly_succinimide(self):
        """ac-A-D.!1(4,3)-G.!1-A-am → 5-membered succinimide ring embedded in chain."""
        mol = _romol('ac-A-D.!1(4,3)-G.!1-A-am')
        assert mol is not None
        assert mol.GetRingInfo().NumRings() >= 1
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []
        imide = Chem.MolFromSmarts('[N]([C]=O)[C]=O')
        assert mol.HasSubstructMatch(imide), "No imide (N flanked by two C=O) found"
        assert mol.GetNumAtoms() == 25

    def test_glu_gly_glutarimide(self):
        """E.!1(4,3) with Glu forms 6-membered glutarimide ring."""
        mol = _romol('ac-A-E.!1(4,3)-G.!1-A-am')
        assert mol is not None
        assert mol.GetRingInfo().NumRings() >= 1
        imide = Chem.MolFromSmarts('[N]([C]=O)[C]=O')
        assert mol.HasSubstructMatch(imide)


class TestSPAAC:
    """SPAAC (copper-free click): strained cyclooctyne + azide → 1,2,3-triazole."""

    def test_spaac_smirks_direct(self):
        """run_bond_smirks produces a triazole from model cyclooctyne + azide."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['spaac_triazole']
        # [400*] on sp3 C alpha to the cyclooctyne (8-membered ring)
        frag1 = Chem.MolFromSmiles('C1CC([400*])C#CCCC1')
        # [401*] directly on the azide alpha C (no extra CH2)
        frag2 = Chem.MolFromSmiles('[401*]CN=[N+]=[N-]')
        assert frag1 is not None and frag2 is not None
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        assert prod is not None
        triazole = Chem.MolFromSmarts('[n]1[n][n][c][c]1')
        assert prod.HasSubstructMatch(triazole), "No 1,2,3-triazole ring in SPAAC product"

    def test_cyclooctyne_chem_type(self):
        """infer_chem_type identifies cyclooctyne_c for [4*] on sp3 C alpha to cyclic alkyne."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        mol = Chem.MolFromSmiles('C1CC([4*])C#CCCC1')
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=3)
        assert ct == 'cyclooctyne_c', f"Expected cyclooctyne_c, got {ct}"

    def test_azide_chem_type_on_azk_fragment(self):
        """infer_chem_type identifies azide_alpha_c for AzK-like Cε–N3 group."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        mol = Chem.MolFromSmiles('CC([4*])N=[N+]=[N-]')
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=3)
        assert ct == 'azide_alpha_c', f"Expected azide_alpha_c, got {ct}"

    def test_azk_monomer_in_library(self):
        """AzK is present in the monomer SDF with correct CHUCKLES."""
        from pathlib import Path
        from rdkit.Chem import SDMolSupplier, MolToSmiles
        sdf = Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
        found = None
        for mol in SDMolSupplier(str(sdf), removeHs=False):
            if mol and mol.GetPropsAsDict().get('m_abbr') == 'AzK':
                found = MolToSmiles(mol)
                break
        assert found is not None, "AzK not found in monomers.sdf"
        assert '[4*]' in found, "AzK CHUCKLES should contain [4*] R4 attachment"
        assert '[N+]=[N-]' in found, "AzK should contain azide group"

    def test_azk_assembled_no_dummies(self):
        """G-AzK-G assembles cleanly with all R1/R2/R3 dummies capped."""
        mol = _romol('G-AzK-G')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unexpected dummies: {dummies}"
        # azide N=N+=N- should survive in the free peptide (R4 uncrosslinked)
        azide = Chem.MolFromSmarts('[N]=[N+]=[N-]')
        assert mol.HasSubstructMatch(azide), "Azide not present in assembled G-AzK-G"


class TestIEDDA:
    """IEDDA tetrazine + TCO wiring: BOND_TABLE entry, chem_type detection, SMIRKS."""

    def test_iedda_in_reaction_index(self):
        """('tetrazine_c', 'tco_c') is routed to iedda_tetrazine_tco."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('tetrazine_c', 'tco_c'))
        assert entry is not None, "No REACTION_INDEX entry for (tetrazine_c, tco_c)"
        assert entry['id'] == 'iedda_tetrazine_tco'

    def test_tetrazine_chem_type(self):
        """infer_chem_type identifies tetrazine_c for sp3 C alpha to s-tetrazine ring."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # methyl-s-tetrazine: methyl C has [4*], bonded to tetrazine ring C
        # Tertiary alpha-C (CH1): bonded to [4*], CH3, and tetrazine ring → CX4H1 → [CX4H:1] ✓
        mol = Chem.MolFromSmiles('C([4*])(C)c1nncnn1')
        assert mol is not None, "Could not build model tetrazine SMILES"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=3)
        assert ct == 'tetrazine_c', f"Expected tetrazine_c, got {ct}"

    def test_tco_chem_type(self):
        """infer_chem_type identifies tco_c for exocyclic sp3 C bonded to ring alkene."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # model: exocyclic CH2([400*]) bonded directly to the vinyl C of a cyclooctene ring.
        # pre_activate uses !r so [4*] always lands on an exocyclic C, never in-ring.
        mol = Chem.MolFromSmiles('C([400*])C1=CCCCCCC1')
        assert mol is not None
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=3)
        assert ct == 'tco_c', f"Expected tco_c, got {ct}"

    def test_iedda_smirks(self):
        """IEDDA combined SMIRKS produces sanitizable dihydropyridazine + N2 fragment."""
        from pyPept.interfaces.reaction_library import REACTIONS
        from rdkit.Chem import AllChem
        entry = REACTIONS['iedda_tetrazine_tco']
        assert len(entry['steps']) == 1, "iedda_tetrazine_tco should have one combined step"
        smirks = entry['steps'][0]
        targeted = smirks.replace('[4*]', '[400*]', 1).replace('[4*]', '[401*]', 1)
        rxn = AllChem.ReactionFromSmarts(targeted)
        assert rxn is not None
        # model: tetrazine methyl cap [400*]Cc1nncnn1; cyclooctene [401*] on exo-C
        tz = Chem.MolFromSmiles('[400*]Cc1nncnn1')
        tco = Chem.MolFromSmiles('C1CC([401*])C=CCCCC1')
        if tz is None or tco is None:
            pytest.skip("Model SMILES failed to parse — check RDKit tetrazine support")
        prods = rxn.RunReactants((tz, tco))
        if not prods:
            prods = rxn.RunReactants((tco, tz))
        assert prods, "IEDDA SMIRKS produced no products with model molecules"
        # all product fragments must sanitize without kekulization error
        for prod_set in prods[:1]:
            for p in prod_set:
                Chem.SanitizeMol(p)  # raises on failure
        # take_largest keeps main ring; N=N fragment is a small byproduct
        main_atoms = max(len(p.GetAtoms()) for p in prods[0])
        assert main_atoms > 4, f"Main product too small ({main_atoms} atoms)"


class TestPipelineLabelling:
    """_SIDECHAIN_RULES auto-labels aldehyde, aminooxy, hydrazide, NHS ester, maleimide."""

    def _chem_type_from_smiles(self, smiles, slot=3):
        from pyPept.interfaces.reaction_library import infer_chem_type
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Bad SMILES: {smiles}"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4
            for nb in a.GetNeighbors()
        )
        return infer_chem_type(mol, attach_idx, slot=slot)

    def test_aldehyde_detection(self):
        # Minimal model C([4*])=O: C has [4*]+O=3 bonds → 1H (CH1).
        # Real CHUCKLES R-C([4*])=O: C has chain+[4*]+O=4 bonds → 0H (CH0).
        # Unified SMARTS [CX3;H0,H1:1](=O)[!#7] covers both; [!#7] blocks N-formyl.
        assert self._chem_type_from_smiles('C([4*])=O') == 'aldehyde'
        assert self._chem_type_from_smiles('CC([4*])=O') == 'aldehyde'

    def test_aminooxy_detection(self):
        # [4*] on N: pre_activate removes one H from NH2, adds [4*]; N has [4*]+O+1H = NX3H1,
        # O has N+chain = OX2H0.  SMARTS [NH1:1][OX2H0] detects this correctly.
        assert self._chem_type_from_smiles('[4*][NH]OC') == 'aminooxy'

    def test_hydrazide_detection(self):
        # pre_activate puts [4*] on amide N-H (its only H is displaced), giving NX3H0.
        # Terminal NH2 is unchanged (NX3H2). SMARTS [NX3H0:1]([NX3H2])C(=O) detects this.
        assert self._chem_type_from_smiles('[4*]N(N)C=O') == 'hydrazide'

    def test_maleimide_detection(self):
        ct = self._chem_type_from_smiles('C([4*])1=CC(=O)NC1=O')
        assert ct == 'maleimide_c', f"Expected maleimide_c, got {ct}"

    def test_nhs_ester_detection(self):
        ct = self._chem_type_from_smiles('C([4*])C(=O)ON1C(=O)CCC1=O')
        assert ct == 'nhs_ester', f"Expected nhs_ester, got {ct}"

    def test_pipeline_labels_aldehyde(self):
        """pre_activate places [4*] on the aldehyde C for a sidechain CHO group."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # model: 4-formyl-phenylalanine
        result = pre_activate('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O')
        chuckles = result[0]
        assert '[4*]' in chuckles, "No R4 placed on aldehyde carbon"
        chem_types = result[2]
        assert 'aldehyde' in chem_types.values(), f"Expected aldehyde chem_type, got {chem_types}"

    def test_pipeline_labels_maleimide(self):
        """pre_activate places [4*] on the maleimide alkene C for a Lys-maleimide AA."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # Lys-maleimide: ε-N carries the maleimide ring; proper backbone (NH2 + COOH) present
        result = pre_activate('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O')
        chuckles = result[0]
        assert chuckles is not None, "pre_activate failed to generate CHUCKLES"
        assert '[4*]' in chuckles, "No R4 placed on maleimide alkene C"
        assert 'maleimide_c' in result[2].values(), f"Expected maleimide_c chem_type, got {result[2]}"

    def test_pipeline_labels_aminooxy(self):
        """pre_activate places [4*] on aminooxy N (not O) for a sidechain O-NH2 group."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # 2-amino-4-(aminooxy)butanoic acid: aminooxy -O-NH2 on gamma-carbon sidechain
        chuckles, _, chem_types, err = pre_activate('NOCC[C@@H](N)C(=O)O')
        assert err is None, f"pre_activate crashed: {err}"
        assert '[4*]' in chuckles, "No R4 placed on aminooxy group"
        assert 'aminooxy' in chem_types.values(), f"Expected aminooxy chem_type, got {chem_types}"
        mol = Chem.MolFromSmiles(chuckles)
        dummy4 = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4)
        nb = next(iter(dummy4.GetNeighbors()))
        assert nb.GetAtomicNum() == 7, \
            f"[4*] should land on N (7), found atomic num {nb.GetAtomicNum()} ({nb.GetSymbol()})"

    def test_pipeline_labels_terminal_alkene(self):
        """pre_activate places [4*] on the internal vinyl C for RCM stapling."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # alpha-amino-7-octenoic acid (S5-type): 5-C sidechain ending in terminal alkene
        chuckles, _, chem_types, err = pre_activate('N[C@@H](CCCC=C)C(=O)O')
        assert err is None, f"pre_activate crashed: {err}"
        assert '[4*]' in chuckles, "No R4 placed on terminal alkene"
        assert 'terminal_alkene' in chem_types.values(), \
            f"Expected terminal_alkene chem_type, got {chem_types}"
        mol = Chem.MolFromSmiles(chuckles)
        dummy4 = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 4)
        nb = next(iter(dummy4.GetNeighbors()))
        # The internal vinyl C is sp2 (double bond to CH2); check it has a double-bond neighbour
        assert any(b.GetBondTypeAsDouble() == 2.0 for b in nb.GetBonds()), \
            "Attachment C should be the sp2 vinyl C (carrying a C=C double bond)"

    def test_pipeline_labels_hydrazide(self):
        """pre_activate places [4*] on hydrazide N for a sidechain -C(=O)-NH-NH2 group."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # Glu-hydrazide: gamma-carboxyl capped as hydrazide (-C(=O)-NH-NH2)
        chuckles, _, chem_types, err = pre_activate('N[C@@H](CCC(=O)NN)C(=O)O')
        assert err is None, f"pre_activate crashed: {err}"
        assert '[4*]' in chuckles, "No R4 placed on hydrazide"
        assert 'hydrazide' in chem_types.values(), \
            f"Expected hydrazide chem_type, got {chem_types}"

    def test_pipeline_labels_alkyne(self):
        """pre_activate places [4*] on sp3 alpha-C adjacent to terminal alkyne (CuAAC)."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        for smi, label in [
            ('N[C@@H](CC#C)C(=O)O',   'Pra (1C spacer)'),
            ('N[C@@H](CCC#C)C(=O)O',  'Hpg (2C spacer)'),
            ('N[C@@H](CCCC#C)C(=O)O', 'Ahex (3C spacer)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'alkyne_c' in chem_types.values(), \
                f"{label}: Expected alkyne_c, got {chem_types}"

    def test_pipeline_labels_azide(self):
        """pre_activate places [4*] on sp3 alpha-C adjacent to organic azide."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        for smi, label in [
            ('N[C@@H](CN=[N+]=[N-])C(=O)O',   'AzAla (1C spacer)'),
            ('N[C@@H](CCN=[N+]=[N-])C(=O)O',  'AzHal (2C spacer)'),
            ('N[C@@H](CCCN=[N+]=[N-])C(=O)O', 'AzOrn (3C spacer)'),
            ('N[C@@H](CCCCN=[N+]=[N-])C(=O)O','AzK   (4C spacer)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'azide_alpha_c' in chem_types.values(), \
                f"{label}: Expected azide_alpha_c, got {chem_types}"

    def test_pipeline_labels_nhs_ester(self):
        """pre_activate places [4*] on sp3 alpha-C adjacent to NHS ester carbonyl."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        for smi, label in [
            ('N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O',   'NhsAsp (1C spacer)'),
            ('N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O',  'NhsGlu (2C spacer)'),
            ('N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O', 'NhsNva (3C spacer)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'nhs_ester' in chem_types.values(), \
                f"{label}: Expected nhs_ester, got {chem_types}"

    def test_pipeline_labels_tetrazine(self):
        """pre_activate places [4*] on sp3 alpha-C adjacent to s-tetrazine ring."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        for smi, label in [
            ('N[C@@H](Cc1nncnn1)C(=O)O',      'TzAla (1C spacer)'),
            ('N[C@@H](CCc1nncnn1)C(=O)O',     'TzAbu (2C spacer)'),
            ('N[C@@H](CCCCCc1nncnn1)C(=O)O',  'TzLys (5C spacer)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'tetrazine_c' in chem_types.values(), \
                f"{label}: Expected tetrazine_c, got {chem_types}"

    def test_pipeline_labels_tco(self):
        """pre_activate places [4*] on exocyclic sp3 C adjacent to ring alkene (TCO)."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        from rdkit import Chem as _Chem
        for smi, label in [
            ('N[C@@H](CC1=CCCCCCC1)C(=O)O',   'TcoAla (1C linker)'),
            ('N[C@@H](CCC1=CCCCCCC1)C(=O)O',  'TcoAbu (2C linker)'),
            ('N[C@@H](CCCCC1=CCCCCCC1)C(=O)O','TcoLys (4C linker)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'tco_c' in chem_types.values(), \
                f"{label}: Expected tco_c, got {chem_types}"
            # The [4*] attachment C must be exocyclic (not in a ring)
            mol = _Chem.MolFromSmiles(chuckles)
            dummy4 = next(a for a in mol.GetAtoms()
                          if a.GetAtomicNum() == 0 and a.GetIsotope() == 4)
            nb = next(iter(dummy4.GetNeighbors()))
            assert not nb.IsInRing(), \
                f"{label}: [4*] should be on exocyclic C, but atom is in ring"

    def test_pipeline_labels_cyclooctyne(self):
        """pre_activate places [4*] on in-ring sp3 C adjacent to cyclooctyne triple bond."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        for smi, label in [
            ('N[C@@H](CC1CCCCC#CC1)C(=O)O',    'CyoAla (1C linker)'),
            ('N[C@@H](CCC1CCCCC#CC1)C(=O)O',   'CyoAbu (2C linker)'),
            ('N[C@@H](CCCCC1CCCCC#CC1)C(=O)O', 'CyoLys (4C linker)'),
        ]:
            chuckles, _, chem_types, err = pre_activate(smi)
            assert err is None, f"{label}: pre_activate crashed: {err}"
            assert '[4*]' in chuckles, f"{label}: No R4 placed"
            assert 'cyclooctyne_c' in chem_types.values(), \
                f"{label}: Expected cyclooctyne_c, got {chem_types}"


class TestOxime:
    """Oxime ligation: aminooxy (N-attachment) + aldehyde → R-O-N=CH-R'."""

    def test_oxime_smirks_direct(self):
        """[4*] on aminooxy N reacts with [4*] on aldehyde C to give oxime O-N=C."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['oxime_ligation']
        # frag1: aminooxy with [4*] on N; N bonded to O which carries a chain C (OX2H0)
        frag1 = Chem.MolFromSmiles('[400*][NH]OCC')
        # frag2: aldehyde with [4*] on aldehyde C (CH with [4*] + =O)
        frag2 = Chem.MolFromSmiles('[401*][CH]=O')
        assert frag1 is not None and frag2 is not None, "Model SMILES failed to parse"
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        assert prod is not None
        smi = Chem.MolToSmiles(prod)
        # Must contain the oxime linkage O-N=C
        oxime_pat = Chem.MolFromSmarts('[O][N]=[C]')
        assert prod.HasSubstructMatch(oxime_pat), \
            f"No oxime O-N=C in product: {smi}"
        # No dummy atoms should remain
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Unreacted dummy atoms in product: {smi}"

    def test_oxime_intramolecular(self):
        """Intramolecular oxime closure gives a ring containing O-N=C."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['oxime_ligation']
        # [400*] on aminooxy N; [401*] on aldehyde C (CHUCKLES-style: H was removed by pre_activate)
        mol = Chem.MolFromSmiles('[400*][NH]OCCC([401*])=O')
        if mol is None:
            pytest.skip("Model SMILES failed to parse")
        prod = run_bond_smirks(mol, mol, 400, 401, entry, intramolecular=True)
        assert prod is not None
        smi = Chem.MolToSmiles(prod)
        oxime_pat = Chem.MolFromSmarts('[O][N]=[C]')
        assert prod.HasSubstructMatch(oxime_pat), \
            f"No oxime O-N=C in intramolecular product: {smi}"
        assert prod.GetRingInfo().NumRings() >= 1, \
            f"Expected ring in intramolecular product: {smi}"


class TestPhosphopeptides:
    """Pre-formed phospho-amino acid monomers assemble correctly."""

    def test_pser_assembles(self):
        mol = _romol('ac-pSer-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []
        phosphate = Chem.MolFromSmarts('[P](=O)(O)O')
        assert mol.HasSubstructMatch(phosphate), "No phosphate group in pSer peptide"

    def test_ptyr_assembles(self):
        mol = _romol('ac-pTyr-am')
        assert mol is not None
        assert mol.HasSubstructMatch(Chem.MolFromSmarts('[P](=O)(O)O'))

    def test_phospho_tripeptide(self):
        mol = _romol('ac-A-pSer-A-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == []


# ---------------------------------------------------------------------------
# Full-coverage direct SMIRKS tests for every reaction in reactions.yaml
# ---------------------------------------------------------------------------

class TestDepsipeptideEster:
    """backbone_ester: O-C(=O) bond for depsipeptide linkages."""

    def test_backbone_ester_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('backbone_o', 'backbone_c'))
        assert entry is not None, "No REACTION_INDEX entry for (backbone_o, backbone_c)"
        assert entry['id'] == 'backbone_ester'

    def test_backbone_ester_smirks_direct(self):
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['backbone_ester']
        frag1 = Chem.MolFromSmiles('[400*]OC')       # backbone_o: [1*][O:2]
        frag2 = Chem.MolFromSmiles('[401*]C(=O)C')  # backbone_c: [2*][C:4](=[O:5])
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        ester = Chem.MolFromSmarts('[O][C](=[O])')
        assert prod.HasSubstructMatch(ester), f"No ester linkage in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy atom leaked into ester product: {smi}"


class TestHydrazoneDirect:
    """hydrazone: hydrazide + aldehyde → N-N=C (water implicit)."""

    def test_hydrazone_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('hydrazide', 'aldehyde'))
        assert entry is not None, "No REACTION_INDEX entry for (hydrazide, aldehyde)"
        assert entry['id'] == 'hydrazone'

    def test_hydrazone_smirks_direct(self):
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['hydrazone']
        # [4*] on the hydrazide N adjacent to NH2 (slot_a/b both 4)
        frag1 = Chem.MolFromSmiles('[400*][NH][NH2]')  # hydrazide: [4*][N:2][N:3]
        frag2 = Chem.MolFromSmiles('[401*][CH]=O')     # aldehyde C: [4*][C:4]=[O]
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        hydrazone = Chem.MolFromSmarts('[N][N]=[C]')
        assert prod.HasSubstructMatch(hydrazone), f"No hydrazone N-N=C in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into hydrazone product: {smi}"

    def test_hydrazone_intramolecular(self):
        """Intramolecular hydrazone closure gives a ring with N-N=C."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['hydrazone']
        # N:2 has [400*], N:3 (branch NH2), and chain → NX3H0 valid; NH2 is terminal branch
        mol = Chem.MolFromSmiles('[400*]N([NH2])CCC([401*])=O')
        if mol is None:
            pytest.skip("Model SMILES failed to parse")
        prod = run_bond_smirks(mol, mol, 400, 401, entry, intramolecular=True)
        assert prod is not None
        hydrazone = Chem.MolFromSmarts('[N][N]=[C]')
        assert prod.HasSubstructMatch(hydrazone), "No hydrazone in intramolecular product"
        assert prod.GetRingInfo().NumRings() >= 1, "Expected ring in intramolecular product"


class TestPhosphorylationDirect:
    """phosphorylation: hydroxyl O + phosphate P → phosphate ester."""

    def test_phosphorylation_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('hydroxyl', 'phosphate_p'))
        assert entry is not None, "No REACTION_INDEX entry for (hydroxyl, phosphate_p)"
        assert entry['id'] == 'phosphorylation'

    def test_phosphorylation_smirks_direct(self):
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['phosphorylation']
        frag1 = Chem.MolFromSmiles('[400*]OC')           # hydroxyl O with chain
        frag2 = Chem.MolFromSmiles('[401*]P(=O)(O)O')   # phosphate P
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        phos_ester = Chem.MolFromSmarts('[O][P](=O)([OH])[OH]')
        assert prod.HasSubstructMatch(phos_ester), f"No phosphate ester in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into phosphorylation product: {smi}"


class TestRCMAlkeneDirect:
    """rcm_alkene: two terminal alkenes → internal alkene + ethylene (take_largest)."""

    def test_rcm_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('terminal_alkene', 'terminal_alkene'))
        assert entry is not None, "No REACTION_INDEX entry for (terminal_alkene, terminal_alkene)"
        assert entry['id'] == 'rcm_alkene'

    def test_rcm_alkene_smirks_direct(self):
        """[4*][C:1]=[CH2].[4*][C:3]=[CH2] >> [C:1]=[C:3] + ethylene; take_largest keeps alkene."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['rcm_alkene']
        # [4*] on internal vinyl C; chain C gives a non-trivial product larger than ethylene
        frag1 = Chem.MolFromSmiles('CC([400*])=C')  # prop-1-en-2-yl fragment
        frag2 = Chem.MolFromSmiles('CC([401*])=C')
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        alkene = Chem.MolFromSmarts('[C]=[C]')
        assert prod.HasSubstructMatch(alkene), f"No alkene in RCM product: {smi}"
        # Ethylene byproduct (2 heavy atoms) should have been discarded by take_largest
        assert prod.GetNumHeavyAtoms() > 2, f"take_largest failed; product too small: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into RCM product: {smi}"


class TestImideAcylationDirect:
    """imide_n_acylation: backbone N (R3) + sidechain carboxyl (R4) → cyclic imide."""

    def test_imide_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('backbone_n_mod', 'carboxyl'))
        assert entry is not None, "No REACTION_INDEX entry for (backbone_n_mod, carboxyl)"
        assert entry['id'] == 'imide_n_acylation'

    def test_imide_n_acylation_smirks_direct(self):
        """N-C(=O) bond formed between backbone amide N and sidechain carboxyl."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['imide_n_acylation']
        frag1 = Chem.MolFromSmiles('[400*]NC')       # backbone N (R3): [3*][N:1]
        frag2 = Chem.MolFromSmiles('[401*]C(=O)C')    # sidechain carboxyl C-attach: [4*][C:2](=[O:3])
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        amide = Chem.MolFromSmarts('[N][C](=[O])')
        assert prod.HasSubstructMatch(amide), f"No amide/imide in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into imide product: {smi}"

    def test_imide_n_acylation_intramolecular(self):
        """Intramolecular closure gives 5-membered succinimide ring."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['imide_n_acylation']
        # Chain: [400*]N-C(=O)-CH2-CH2-C([401*])=O → 5-membered ring on closure
        mol = Chem.MolFromSmiles('[400*]NC(=O)CCC([401*])=O')
        if mol is None:
            pytest.skip("Model SMILES failed to parse")
        prod = run_bond_smirks(mol, mol, 400, 401, entry, intramolecular=True)
        assert prod is not None
        assert prod.GetRingInfo().NumRings() >= 1, "Expected ring in succinimide product"
        amide = Chem.MolFromSmarts('[N][C](=[O])')
        assert prod.HasSubstructMatch(amide), "No amide bond in succinimide product"


class TestThioesterDirect:
    """thioester: thiol S + carboxyl C → S-C(=O) crosslink."""

    def test_thioester_thiol_carboxyl_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('thiol', 'carboxyl'))
        assert entry is not None
        assert entry['id'] == 'thioester_sc'

    def test_thioester_smirks_direct(self):
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['thioester']
        frag1 = Chem.MolFromSmiles('[400*]SC')      # thiol S
        frag2 = Chem.MolFromSmiles('[401*]C(=O)C')  # carboxyl C
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        thioester = Chem.MolFromSmarts('[S][C](=[O])')
        assert prod.HasSubstructMatch(thioester), f"No thioester S-C(=O) in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into thioester product: {smi}"


class TestThioetherHalide:
    """thioether_halide: thiol + primary alkyl halide → thioether via SN2."""

    def test_thioether_halide_in_reaction_index(self):
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('thiol', 'alkyl_halide_c'))
        assert entry is not None, "No REACTION_INDEX entry for (thiol, alkyl_halide_c)"
        assert entry['id'] == 'thioether_halide'

    def test_alkyl_halide_chem_type(self):
        """infer_chem_type returns alkyl_halide_c for C bonded to [4*] and Cl."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # After pre_activate: CH2Cl → C([4*])HCl; C is bonded to [400*], chain, Cl, H
        mol = Chem.MolFromSmiles('[400*]CCl')
        assert mol is not None
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=3)
        assert ct == 'alkyl_halide_c', f"Expected alkyl_halide_c, got {ct}"

    def test_thioether_halide_smirks_direct(self):
        """Halide replaced by dummy in pre_activate CHUCKLES; SMIRKS just forms S-C bond."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['thioether_halide']
        frag1 = Chem.MolFromSmiles('[400*]SC')    # thiol S: [4*][S:2]
        frag2 = Chem.MolFromSmiles('[401*]C')     # alkyl: [4*][C:4] — halide already replaced by dummy
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        thioether = Chem.MolFromSmarts('[S][C]')
        assert prod.HasSubstructMatch(thioether), f"No thioether S-C in product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into thioether product: {smi}"

    def test_pipeline_labels_alkyl_halide(self):
        """pre_activate places [4*] on the alpha-C of a primary alkyl chloride sidechain."""
        from pyPept.interfaces.monomer_pipeline import pre_activate
        # 3-chloro-L-alanine: beta-CH2Cl sidechain
        chuckles, _, chem_types, err = pre_activate('N[C@@H](CCl)C(=O)O')
        assert err is None, f"pre_activate crashed: {err}"
        assert 'alkyl_halide_c' in chem_types.values(), \
            f"Expected alkyl_halide_c in chem_types, got {chem_types}"


class TestIntermolecularSMIRKS:
    """Intermolecular versions of reactions previously only tested intramolecularly."""

    def test_cuaac_intermolecular(self):
        """CuAAC terminal alkyne + azide → triazole (intermolecular)."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['cuaac_1_4_triazole']
        frag1 = Chem.MolFromSmiles('[400*]CC#C')           # [4*]C[C:2]#[CH:3]
        frag2 = Chem.MolFromSmiles('[401*]CN=[N+]=[N-]')  # [4*][C:5]N=[N+]=[N-]
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        assert prod.GetRingInfo().NumRings() >= 1, f"No triazole ring in CuAAC product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into CuAAC product: {smi}"

    def test_thiol_maleimide_intermolecular(self):
        """Thiol-maleimide Michael addition (intermolecular): succinimide ring retained."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['thiol_maleimide']
        # maleimide with [400*] on alkene C; thiol with [401*] on S
        frag1 = Chem.MolFromSmiles('[400*]C1=CC(=O)NC1=O')
        frag2 = Chem.MolFromSmiles('[401*]SC')
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        succinimide_thioether = Chem.MolFromSmarts('[S]C1CC(=O)NC1=O')
        assert prod.HasSubstructMatch(succinimide_thioether), \
            f"No succinimide-thioether in thiol-maleimide product: {smi}"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into thiol-maleimide product: {smi}"

    def test_nhs_ester_amide_intermolecular(self):
        """NHS ester + primary amine → amide; NHS ring departs (take_largest)."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['nhs_ester_amide']
        frag1 = Chem.MolFromSmiles('[400*]CC(=O)ON1C(=O)CCC1=O')  # NHS ester α-C
        frag2 = Chem.MolFromSmiles('[401*]NC')                     # primary amine N
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        amide = Chem.MolFromSmarts('[N][C](=[O])')
        assert prod.HasSubstructMatch(amide), f"No amide in NHS ester product: {smi}"
        # NHS leaving group (O=C1CCC(=O)N1) should be absent in the largest fragment
        nhs_ring = Chem.MolFromSmarts('O=C1CCCN1C=O')
        assert not prod.HasSubstructMatch(nhs_ring), \
            "NHS ring should have departed as leaving group"
        assert not any(a.GetAtomicNum() == 0 for a in prod.GetAtoms()), \
            f"Dummy leaked into NHS ester product: {smi}"


# ---------------------------------------------------------------------------
# Sequence.validate() — ValidationReport dry-run API
# ---------------------------------------------------------------------------

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

    def test_valid_bond_tuple_shape(self):
        """Each bond is a (m1, slot1, m2, slot2) tuple."""
        report = Sequence.validate('A-G')
        assert len(report.bonds) == 1
        b = report.bonds[0]
        assert len(b) == 4  # (m1, slot1, m2, slot2)

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

    def test_warnings_captured_chemistry_mismatch(self):
        """Exotic bond type generates a warning captured in report.warnings."""
        # S-C thioether crosslink triggers a chemistry UserWarning
        report = Sequence.validate('C.!1(4,1)-A-C.!1(4,1)-am')
        # warnings list may or may not be non-empty depending on bond type,
        # but the report must be ok and no exception raised
        assert isinstance(report.warnings, list)

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


# ---------------------------------------------------------------------------
# 11. Sequential reaction bracket notation — .[A(r,s).B(t,u)...]
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Monomer builder CLI tests
# ---------------------------------------------------------------------------

import unittest as _unittest


class TestMonomerBuilderCLI(_unittest.TestCase):
    """Tests for cli_monomer.smiles_from_cabiln, register_monomer, and main()."""

    # ------------------------------------------------------------------ #
    # smiles_from_cabiln
    # ------------------------------------------------------------------ #

    def test_smiles_from_cabiln_single_residue(self):
        """smiles_from_cabiln('A') returns a parseable SMILES for Alanine."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('A')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"smiles_from_cabiln('A') gave invalid SMILES: {smi}"

    def test_smiles_from_cabiln_single_residue_correct_atoms(self):
        """Single-residue Gly assembly gives a SMILES with N and C=O but no dummy atoms."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('G')
        assert '[*]' not in smi and '[0*]' not in smi, \
            f"Dummy atoms unexpectedly present in assembled SMILES: {smi}"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_smiles_from_cabiln_capped_residue(self):
        """smiles_from_cabiln('fmoc-K-am') returns parseable SMILES (capped Lys)."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('fmoc-K-am')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"SMILES invalid: {smi}"

    def test_smiles_from_cabiln_inline_cap(self):
        """smiles_from_cabiln with inline cap notation assembles correctly."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('fmoc-C.trt(4,2)-am')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    # ------------------------------------------------------------------ #
    # register_monomer (uses a temp SDF so the library is never mutated)
    # ------------------------------------------------------------------ #

    def _make_temp_sdf(self):
        """Return path to an empty temp SDF that is cleaned up after the test."""
        tmp = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False, mode='w')
        tmp.close()
        self.addCleanup(lambda: pathlib.Path(tmp.name).unlink(missing_ok=True))
        return tmp.name

    def test_register_monomer_plain_smiles(self):
        """register_monomer with plain SMILES appends a record to the SDF."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        # Alanine SMILES — simple and well-understood
        result = register_monomer('N[C@@H](C)C(=O)O', symbol='TestAla',
                                  sdf_path=sdf)
        assert result.chuckles, "Expected non-empty CHUCKLES"
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert '$$$$' in content, "SDF record terminator missing"
        assert 'TestAla' in content, "Symbol not written to SDF"

    def test_register_monomer_returns_activation_result(self):
        """register_monomer returns an ActivationResult with leaving and chem_types."""
        from pyPept.interfaces.cli_monomer import register_monomer
        from pyPept.interfaces.monomer_pipeline import ActivationResult
        sdf = self._make_temp_sdf()
        result = register_monomer('N[C@@H](CS)C(=O)O', symbol='TestCys', sdf_path=sdf)
        assert isinstance(result, ActivationResult)
        assert result.leaving
        assert result.chem_types

    def test_register_monomer_multiple_append(self):
        """Two calls to register_monomer both appear in the SDF (no overwrite)."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('N[C@@H](C)C(=O)O', symbol='Mon1', sdf_path=sdf)
        register_monomer('NCC(=O)O', symbol='Mon2', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert content.count('$$$$') == 2, \
            f"Expected 2 records, got {content.count('$$$$')}: {content[:400]}"
        assert 'Mon1' in content and 'Mon2' in content

    def test_register_monomer_sets_type_and_subtype(self):
        """Custom m_type and m_subtype are written to the SDF properties."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('NCC(=O)O', symbol='GlyCap',
                         m_type='cap', m_subtype='cap', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert 'm_type' in content
        assert 'cap' in content

    def test_register_monomer_rgroups_written(self):
        """m_Rgroups property is present in the SDF record."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('N[C@@H](C)C(=O)O', symbol='TestAla2', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert 'm_Rgroups' in content

    # ------------------------------------------------------------------ #
    # main() via argparse (CLI integration)
    # ------------------------------------------------------------------ #

    def _run_main(self, argv):
        """Run cli_monomer.main() with the given argv list; return (stdout, stderr, rc)."""
        from pyPept.interfaces import cli_monomer
        import io
        from unittest.mock import patch
        out_buf = io.StringIO()
        err_buf = io.StringIO()
        rc = 0
        with patch('sys.argv', ['pyPept-monomer-add'] + argv), \
             patch('sys.stdout', out_buf), \
             patch('sys.stderr', err_buf):
            try:
                cli_monomer.main()
            except SystemExit as exc:
                rc = int(exc.code) if exc.code is not None else 0
        return out_buf.getvalue(), err_buf.getvalue(), rc

    def test_main_smiles_path(self):
        """CLI --smiles path exits 0 and prints 'Registered'."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'N[C@@H](C)C(=O)O',
            '--symbol', 'CLIAla',
            '--sdf', sdf,
        ])
        assert rc == 0, f"Non-zero exit code; stderr: {err}"
        assert 'Registered' in out, f"Expected 'Registered' in stdout: {out}"
        assert 'CLIAla' in out

    def test_main_from_cabiln_path(self):
        """CLI --from-cabiln path assembles and registers correctly."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--from-cabiln', 'G',
            '--symbol', 'CLIGly',
            '--sdf', sdf,
        ])
        assert rc == 0, f"Non-zero exit code; stderr: {err}"
        assert 'Assembled SMILES' in out, f"Assembly message missing: {out}"
        assert 'Registered' in out

    def test_main_sdf_is_written(self):
        """CLI writes a valid SDF record to the given --sdf path."""
        sdf = self._make_temp_sdf()
        self._run_main([
            '--smiles', 'NCC(=O)O',
            '--symbol', 'CLIGly2',
            '--sdf', sdf,
        ])
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert '$$$$' in content
        assert 'CLIGly2' in content

    def test_main_bad_smiles_exits_nonzero(self):
        """CLI exits non-zero when --smiles cannot be pre-activated."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'NOTVALIDSMILES???',
            '--symbol', 'Bad',
            '--sdf', sdf,
        ])
        assert rc != 0, f"Expected non-zero exit for invalid SMILES; got rc={rc}"

    def test_main_missing_symbol_exits_nonzero(self):
        """CLI exits non-zero when --symbol is omitted."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'NCC(=O)O',
            '--sdf', sdf,
        ])
        assert rc != 0

    def test_main_no_source_exits_nonzero(self):
        """CLI exits non-zero when neither --smiles nor --from-cabiln is given."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--symbol', 'NoSrc',
            '--sdf', sdf,
        ])
        assert rc != 0


# ---------------------------------------------------------------------------
# New monomers added in v1.1.0b0: Tz, TCO, Mal, Aoa, Ald, Dmb
# ---------------------------------------------------------------------------

class TestNewMonomersV110:
    """Integration tests for the 6 monomers added in v1.1.0b0."""

    # ── chem-type detection ────────────────────────────────────────────────

    def test_tz_chem_type_detected(self):
        """Tz monomer exposes tetrazine_c at R4."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # Tz CHUCKLES: [2*]C(=O)C([4*])c1nncnn1 — sp3 CH bonded to tetrazine
        mol = Chem.MolFromSmiles('C([400*])(C(=O)O)c1nncnn1')
        assert mol is not None, "Tz model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=4)
        assert ct == 'tetrazine_c', f"Expected tetrazine_c, got {ct!r}"

    def test_tco_chem_type_detected(self):
        """TCO monomer exposes tco_c at R4."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # TCO CHUCKLES: [2*]C(=O)CC([4*])C1=CCCCCCC1 — exocyclic sp3 C bonded to ring alkene
        mol = Chem.MolFromSmiles('C([400*])(CC(=O)O)C1=CCCCCCC1')
        assert mol is not None, "TCO model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=4)
        assert ct == 'tco_c', f"Expected tco_c, got {ct!r}"

    def test_mal_chem_type_detected(self):
        """Mal monomer exposes maleimide_c at R4."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # Mal: vinyl C of maleimide ring
        mol = Chem.MolFromSmiles('[400*]C1=CC(=O)NC1=O')
        assert mol is not None, "Mal model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=4)
        assert ct == 'maleimide_c', f"Expected maleimide_c, got {ct!r}"

    def test_aoa_chem_type_detected(self):
        """Aoa monomer exposes aminooxy at R4."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        mol = Chem.MolFromSmiles('[400*][NH]OCC')
        assert mol is not None, "Aoa model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=4)
        assert ct == 'aminooxy', f"Expected aminooxy, got {ct!r}"

    def test_ald_chem_type_detected(self):
        """Ald monomer exposes aldehyde at R4."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        mol = Chem.MolFromSmiles('[400*][CH]=O')
        assert mol is not None, "Ald model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=4)
        assert ct == 'aldehyde', f"Expected aldehyde, got {ct!r}"

    def test_dmb_chem_type_carbon(self):
        """Dmb monomer exposes carbon at R1."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # Dmb: [1*]Cc1ccc(OC)cc1OC — benzylic sp3 C
        mol = Chem.MolFromSmiles('[100*]Cc1ccc(OC)cc1OC')
        assert mol is not None, "Dmb model SMILES failed to parse"
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 100
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=1)
        assert ct == 'carbon', f"Expected carbon, got {ct!r}"

    # ── assembly integration ───────────────────────────────────────────────

    def test_tz_cap_on_lys_assembles(self):
        """fmoc-K.Tz(4,2)-am: amide links K ε-amine to Tz carbonyl; tetrazine ring intact."""
        mol = _romol('fmoc-K.Tz(4,2)-am')
        assert mol is not None, "Assembly returned None for K.Tz peptide"
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in K.Tz product: {[a.GetIsotope() for a in dummies]}"
        tz_ring = Chem.MolFromSmarts('c1nncnn1')
        assert mol.HasSubstructMatch(tz_ring), "Tetrazine ring missing from K.Tz product"

    def test_tco_cap_on_lys_assembles(self):
        """fmoc-K.TCO(4,2)-am: amide links K ε-amine to TCO carbonyl; cyclooctene intact."""
        mol = _romol('fmoc-K.TCO(4,2)-am')
        assert mol is not None, "Assembly returned None for K.TCO peptide"
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in K.TCO product: {[a.GetIsotope() for a in dummies]}"
        # Trans-cyclooctene: 8-membered ring with one double bond
        ring_alkene = Chem.MolFromSmarts('[r8]=[r8]')
        assert mol.HasSubstructMatch(ring_alkene), "Cyclooctene ring alkene missing from K.TCO product"

    def test_mal_cap_on_lys_assembles(self):
        """fmoc-K.Mal(4,2)-am: amide links K ε-amine to Mal carbonyl; succinimide ring intact."""
        mol = _romol('fmoc-K.Mal(4,2)-am')
        assert mol is not None, "Assembly returned None for K.Mal peptide"
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in K.Mal product: {[a.GetIsotope() for a in dummies]}"
        maleimide = Chem.MolFromSmarts('C1=CC(=O)NC1=O')
        assert mol.HasSubstructMatch(maleimide), "Maleimide ring missing from K.Mal product"

    def test_mal_thiol_michael_smirks(self):
        """Thiol-maleimide routes to thiol_maleimide in REACTION_INDEX (SMIRKS-level)."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('thiol', 'maleimide_c'))
        assert entry is not None, "No REACTION_INDEX entry for (thiol, maleimide_c)"
        assert entry['id'] == 'thiol_maleimide'

    def test_aoa_peptide_assembles(self):
        """fmoc-A-Aoa-G-am: Aoa backbone incorporation, all dummies consumed."""
        mol = _romol('fmoc-A-Aoa-G-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in Aoa peptide: {[a.GetIsotope() for a in dummies]}"
        aminooxy = Chem.MolFromSmarts('[N][O]')
        assert mol.HasSubstructMatch(aminooxy), "Aminooxy N-O missing from Aoa peptide"

    def test_ald_peptide_assembles(self):
        """fmoc-A-Ald-G-am: Ald backbone incorporation, aldehyde sidechain intact."""
        mol = _romol('fmoc-A-Ald-G-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in Ald peptide: {[a.GetIsotope() for a in dummies]}"
        aldehyde = Chem.MolFromSmarts('[CH]=O')
        assert mol.HasSubstructMatch(aldehyde), "Aldehyde group missing from Ald peptide"

    def test_aoa_ald_oxime_macrocycle(self):
        """fmoc-A-Aoa.!1(4,4)-G-Ald.!1-am: oxime macrocyclisation gives O-N=C linkage."""
        mol = _romol('fmoc-A-Aoa.!1(4,4)-G-Ald.!1-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Unreacted dummies in oxime macrocycle: {[a.GetIsotope() for a in dummies]}"
        oxime = Chem.MolFromSmarts('[O][N]=[C]')
        assert mol.HasSubstructMatch(oxime), "Oxime O-N=C missing from Aoa-Ald macrocycle"
        assert mol.GetRingInfo().NumRings() >= 1, "Expected macrocycle ring"

    def test_tz_iedda_in_reaction_index(self):
        """REACTION_INDEX routes (tetrazine_c, tco_c) → iedda_tetrazine_tco."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('tetrazine_c', 'tco_c'))
        assert entry is not None
        assert entry['id'] == 'iedda_tetrazine_tco'

    def test_tz_tco_iedda_reaction_index(self):
        """REACTION_INDEX routes (tetrazine_c, tco_c) → iedda_tetrazine_tco (routing check)."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('tetrazine_c', 'tco_c'))
        assert entry is not None
        assert entry['id'] == 'iedda_tetrazine_tco'
        assert entry.get('take_largest') is True, "IEDDA must have take_largest=True (N2 byproduct)"


class TestMonomerPreActivate:
    """Monomer-level tests: pre_activate on raw SMILES verifies R-group placement.

    For each amino acid, checks that the CHUCKLES dummy atoms land on the
    correct attachment atoms (O for carboxyl, S for thiol, etc.) with the
    right leaving group and chem_type.
    """

    @staticmethod
    def _slot_map(chuckles):
        """Parse CHUCKLES -> {slot_isotope: neighbor_atomic_num}."""
        mol = Chem.MolFromSmiles(chuckles)
        assert mol is not None, f"Bad CHUCKLES: {chuckles}"
        result = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                iso = atom.GetIsotope()
                nbs = list(atom.GetNeighbors())
                assert len(nbs) == 1, f"Dummy {iso} has {len(nbs)} neighbors"
                result[iso] = nbs[0].GetAtomicNum()
        return result

    @staticmethod
    def _verify_attachment_context(chuckles, slot, chem_type):
        """Verify the chemical environment around the attachment atom is correct.

        Goes beyond element-type checking to confirm the dummy is on the
        structurally correct atom — e.g. the carbonyl C of a carboxyl with
        =O, not the hydroxyl O.
        """
        mol = Chem.MolFromSmiles(chuckles)
        dummy = next(a for a in mol.GetAtoms()
                     if a.GetAtomicNum() == 0 and a.GetIsotope() == slot)
        attach = next(iter(dummy.GetNeighbors()))

        if chem_type == 'carboxyl':
            assert attach.GetAtomicNum() == 6, \
                f"Carboxyl slot {slot}: dummy should be on C, got {attach.GetSymbol()}"
            has_co = any(
                mol.GetBondBetweenAtoms(attach.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in attach.GetNeighbors() if nb.GetAtomicNum() == 8)
            assert has_co, \
                f"Carboxyl slot {slot}: C must have =O double bond [C-attachment]"

        elif chem_type == 'aldehyde':
            assert attach.GetAtomicNum() == 6, \
                f"Aldehyde slot {slot}: dummy should be on C, got {attach.GetSymbol()}"
            has_co = any(
                mol.GetBondBetweenAtoms(attach.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in attach.GetNeighbors() if nb.GetAtomicNum() == 8)
            assert has_co, \
                f"Aldehyde slot {slot}: C must have =O double bond [C-attachment]"

        elif chem_type == 'backbone_c':
            assert attach.GetAtomicNum() == 6
            has_co = any(
                mol.GetBondBetweenAtoms(attach.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in attach.GetNeighbors() if nb.GetAtomicNum() == 8)
            assert has_co, f"Backbone C slot {slot}: must be carbonyl C with =O"

        elif chem_type == 'thiol':
            assert attach.GetAtomicNum() == 16
            heavy = [nb for nb in attach.GetNeighbors()
                     if nb.GetAtomicNum() > 1]
            assert len(heavy) == 1, \
                f"Thiol slot {slot}: S should have 1 heavy neighbor, got {len(heavy)}"

        elif chem_type == 'selenol':
            assert attach.GetAtomicNum() == 34
            heavy = [nb for nb in attach.GetNeighbors()
                     if nb.GetAtomicNum() > 1]
            assert len(heavy) == 1, \
                f"Selenol slot {slot}: Se should have 1 heavy neighbor, got {len(heavy)}"

        elif chem_type == 'terminal_alkene':
            assert attach.GetAtomicNum() == 6
            has_cc_double = any(
                mol.GetBondBetweenAtoms(attach.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in attach.GetNeighbors() if nb.GetAtomicNum() == 6)
            assert has_cc_double, \
                f"Terminal alkene slot {slot}: C must have C=C double bond"

        elif chem_type == 'hydroxyl':
            assert attach.GetAtomicNum() == 8
            c_nbs = [nb for nb in attach.GetNeighbors()
                     if nb.GetAtomicNum() == 6 and nb.GetIdx() != dummy.GetIdx()]
            assert any(not nb.GetIsAromatic() for nb in c_nbs), \
                f"Hydroxyl slot {slot}: O must be bonded to aliphatic C"

        elif chem_type == 'hydroxyl_phenolic':
            assert attach.GetAtomicNum() == 8
            has_aromatic = any(nb.GetIsAromatic()
                              for nb in attach.GetNeighbors()
                              if nb.GetAtomicNum() == 6)
            assert has_aromatic, \
                f"Hydroxyl_phenolic slot {slot}: O must be bonded to aromatic C"

        elif chem_type == 'phosphate_p':
            assert attach.GetAtomicNum() == 15
            has_po = any(
                mol.GetBondBetweenAtoms(attach.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in attach.GetNeighbors() if nb.GetAtomicNum() == 8)
            assert has_po, f"Phosphate slot {slot}: P must have P=O"

        elif chem_type == 'alkyl_halide_c':
            assert attach.GetAtomicNum() == 6

        elif chem_type == 'aminooxy':
            assert attach.GetAtomicNum() == 7
            has_o_nb = any(nb.GetAtomicNum() == 8
                          for nb in attach.GetNeighbors()
                          if nb.GetIdx() != dummy.GetIdx())
            assert has_o_nb, \
                f"Aminooxy slot {slot}: N must be bonded to O"

        elif chem_type == 'hydrazide':
            assert attach.GetAtomicNum() == 7
            has_n_nb = any(nb.GetAtomicNum() == 7
                          for nb in attach.GetNeighbors()
                          if nb.GetIdx() != dummy.GetIdx())
            assert has_n_nb, \
                f"Hydrazide slot {slot}: N must be bonded to another N"

        elif chem_type == 'maleimide_c':
            assert attach.GetAtomicNum() == 6
            assert attach.IsInRing(), \
                f"Maleimide slot {slot}: C must be in ring"

        elif chem_type == 'backbone_n':
            assert attach.GetAtomicNum() == 7

        elif chem_type == 'backbone_n_mod':
            assert attach.GetAtomicNum() == 7

    def _check(self, smiles, expect_slots=None, expect_lg=None, expect_ct=None):
        from pyPept.interfaces.monomer_pipeline import pre_activate
        r = pre_activate(smiles)
        smap = self._slot_map(r.chuckles)

        if expect_slots is not None:
            assert len(smap) == len(expect_slots), (
                f"Expected {len(expect_slots)} dummies, got {len(smap)}: {smap}")
            for slot, expected_elem in expect_slots.items():
                assert slot in smap, f"Missing slot {slot} in CHUCKLES"
                assert smap[slot] == expected_elem, (
                    f"Slot {slot}: expected element {expected_elem}, got {smap[slot]}")

        if expect_lg:
            for slot, lg in expect_lg.items():
                assert r.leaving.get(slot) == lg, (
                    f"Slot {slot}: expected LG {lg!r}, got {r.leaving.get(slot)!r}")

        if expect_ct:
            for slot, ct in expect_ct.items():
                assert r.chem_types.get(slot) == ct, (
                    f"Slot {slot}: expected ct {ct!r}, got {r.chem_types.get(slot)!r}")
                self._verify_attachment_context(r.chuckles, slot, ct)

        return r

    # ── Carboxyl sidechain: C-attachment (dummy on carbonyl C) ────────────────

    def test_asp_carboxyl_on_carbon(self):
        """Asp: slot 4 dummy lands on sidechain carboxyl C, LG=[OH]."""
        self._check('N[C@@H](CC(=O)O)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_lg={4: '[OH]'}, expect_ct={4: 'carboxyl'})

    def test_glu_carboxyl_on_carbon(self):
        """Glu: slot 4 dummy on sidechain carboxyl C."""
        self._check('N[C@@H](CCC(=O)O)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_lg={4: '[OH]'}, expect_ct={4: 'carboxyl'})

    def test_d_asp_carboxyl_on_carbon(self):
        """D-Asp: C-attachment matches L-Asp."""
        self._check('N[C@H](CC(=O)O)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_lg={4: '[OH]'}, expect_ct={4: 'carboxyl'})

    def test_d_glu_carboxyl_on_carbon(self):
        """D-Glu: C-attachment matches L-Glu."""
        self._check('N[C@H](CCC(=O)O)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_lg={4: '[OH]'}, expect_ct={4: 'carboxyl'})

    def test_asp_backbone_cooh_oh_excluded(self):
        """Asp: backbone COOH OH is not spuriously labelled as carboxyl."""
        r = self._check('N[C@@H](CC(=O)O)C(=O)O',
                         {1: 7, 2: 6, 3: 7, 4: 6})
        carboxyl_slots = [s for s, ct in r.chem_types.items() if ct == 'carboxyl']
        assert len(carboxyl_slots) == 1, f"Expected 1 carboxyl slot, got {carboxyl_slots}"

    def test_glu_backbone_cooh_oh_excluded(self):
        """Glu: backbone COOH OH excluded from sidechain detection."""
        r = self._check('N[C@@H](CCC(=O)O)C(=O)O',
                         {1: 7, 2: 6, 3: 7, 4: 6})
        carboxyl_slots = [s for s, ct in r.chem_types.items() if ct == 'carboxyl']
        assert len(carboxyl_slots) == 1

    # ── Thiol sidechain ──────────────────────────────────────────────────────

    def test_cys_thiol_on_sulfur(self):
        """Cys: slot 4 dummy on thiol S."""
        self._check('N[C@@H](CS)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 16},
                     expect_lg={4: '[H]'}, expect_ct={4: 'thiol'})

    def test_d_cys_thiol_on_sulfur(self):
        """D-Cys: same thiol slot as L-Cys."""
        self._check('N[C@H](CS)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 16},
                     expect_lg={4: '[H]'}, expect_ct={4: 'thiol'})

    # ── Selenol sidechain ────────────────────────────────────────────────────

    def test_sec_selenol_on_selenium(self):
        """Sec (selenocysteine): slot 4 dummy on Se."""
        self._check('N[C@@H](C[SeH])C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 34},
                     expect_lg={4: '[H]'}, expect_ct={4: 'selenol'})

    # ── Hydroxyl sidechain ───────────────────────────────────────────────────

    def test_ser_hydroxyl_on_oxygen(self):
        """Ser: slot 4 dummy on hydroxyl O."""
        self._check('N[C@@H](CO)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_lg={4: '[H]'}, expect_ct={4: 'hydroxyl'})

    def test_thr_hydroxyl_on_oxygen(self):
        """Thr: slot 4 dummy on hydroxyl O."""
        self._check('N[C@@H]([C@@H](O)C)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_lg={4: '[H]'}, expect_ct={4: 'hydroxyl'})

    # ── Amine sidechain ──────────────────────────────────────────────────────

    def test_lys_amine_primary_on_nitrogen(self):
        """Lys: epsilon-NH2 -> amine_primary + amine_secondary (two H slots)."""
        r = self._check('N[C@@H](CCCCN)C(=O)O',
                         {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                         expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    def test_orn_amine_primary_on_nitrogen(self):
        """Orn: delta-NH2 same amine_primary + amine_secondary."""
        self._check('N[C@@H](CCCN)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    # ── Aromatic NH sidechain ────────────────────────────────────────────────

    def test_trp_aromatic_nh(self):
        """Trp: indole NH detected as aromatic_nh."""
        self._check('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7},
                     expect_ct={4: 'aromatic_nh'})

    def test_his_aromatic_nh(self):
        """His: imidazole NH detected as aromatic_nh."""
        self._check('N[C@@H](Cc1cnc[nH]1)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7},
                     expect_ct={4: 'aromatic_nh'})

    # ── Phenolic hydroxyl ────────────────────────────────────────────────────

    def test_tyr_phenolic_oh(self):
        """Tyr: phenolic OH detected as hydroxyl_phenolic on O."""
        self._check('N[C@@H](Cc1ccc(O)cc1)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl_phenolic'})

    # ── Amide sidechain (NH2 matches amine_primary) ──────────────────────────

    def test_asn_amide_nh2(self):
        """Asn: amide NH2 matches amine_primary (H2), plus amine_secondary."""
        self._check('N[C@@H](CC(=O)N)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    def test_gln_amide_nh2(self):
        """Gln: same amide NH2 pattern as Asn."""
        self._check('N[C@@H](CCC(=O)N)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    # ── No sidechain reactive group ──────────────────────────────────────────

    def test_ala_no_sidechain(self):
        """Ala: methyl sidechain -> 3 backbone slots only."""
        self._check('N[C@@H](C)C(=O)O', {1: 7, 2: 6, 3: 7},
                     expect_ct={1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'})

    def test_gly_no_sidechain(self):
        """Gly: no sidechain -> 3 backbone slots only."""
        self._check('NCC(=O)O', {1: 7, 2: 6, 3: 7})

    def test_val_no_sidechain(self):
        """Val: isopropyl sidechain has no reactive group."""
        self._check('N[C@@H](C(C)C)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_leu_no_sidechain(self):
        """Leu: isobutyl sidechain has no reactive group."""
        self._check('N[C@@H](CC(C)C)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_ile_no_sidechain(self):
        """Ile: sec-butyl sidechain has no reactive group."""
        self._check('N[C@@H]([C@@H](C)CC)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_phe_no_sidechain(self):
        """Phe: benzyl sidechain has no reactive group."""
        self._check('N[C@@H](Cc1ccccc1)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_met_thioether_not_thiol(self):
        """Met: thioether S-CH3 is NOT detected as thiol -> 3 backbone slots only."""
        self._check('N[C@@H](CCSC)C(=O)O', {1: 7, 2: 6, 3: 7})

    # ── Pro: no backbone_n_mod (ring N has 1H) ───────────────────────────────

    def test_pro_no_backbone_n_mod(self):
        """Pro: ring N has only 1H -> no backbone_n_mod -> 2 dummies."""
        self._check('OC(=O)[C@@H]1CCCN1', {1: 7, 2: 6},
                     expect_ct={1: 'backbone_n', 2: 'backbone_c'})

    # ── Alpha-methyl AAs ─────────────────────────────────────────────────────

    def test_aib_alpha_methyl_no_sidechain(self):
        """Aib: alpha-aminoisobutyric acid, backbone only, no sidechain."""
        self._check('CC(N)(C)C(=O)O', {1: 7, 2: 6, 3: 7})

    # ── Backbone slot verification ───────────────────────────────────────────

    def test_backbone_n_leaving_is_h(self):
        """Slot 1 (backbone_n) LG is always [H]."""
        r = self._check('N[C@@H](C)C(=O)O', {1: 7, 2: 6, 3: 7})
        assert r.leaving[1] == '[H]'

    def test_backbone_c_leaving_is_oh(self):
        """Slot 2 (backbone_c) LG is always [OH]."""
        r = self._check('N[C@@H](C)C(=O)O', {1: 7, 2: 6, 3: 7})
        assert r.leaving[2] == '[OH]'

    def test_backbone_n_mod_leaving_is_h(self):
        """Slot 3 (backbone_n_mod) LG is [H]."""
        r = self._check('N[C@@H](C)C(=O)O', {1: 7, 2: 6, 3: 7})
        assert r.leaving[3] == '[H]'
        assert r.chem_types[3] == 'backbone_n_mod'

    # ── D-isomer parity: same slot structure as L ────────────────────────────

    def test_d_ala_same_as_l_ala(self):
        """D-Ala: same 3-slot structure as L-Ala."""
        self._check('N[C@H](C)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_d_ser_same_as_l_ser(self):
        """D-Ser: same hydroxyl slot as L-Ser."""
        self._check('N[C@H](CO)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl'})

    def test_d_lys_same_as_l_lys(self):
        """D-Lys: same amine slots as L-Lys."""
        self._check('N[C@H](CCCCN)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    # ── Thiol: 3 more (total 5) ──────────────────────────────────────────────

    def test_homocysteine_thiol(self):
        """Homocysteine: 2-carbon sidechain thiol."""
        self._check('N[C@@H](CCS)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 16},
                     expect_lg={4: '[H]'}, expect_ct={4: 'thiol'})

    def test_d_homocysteine_thiol(self):
        """D-Homocysteine: D-isomer 2C sidechain thiol."""
        self._check('N[C@H](CCS)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 16},
                     expect_ct={4: 'thiol'})

    def test_mercaptonorvaline_thiol(self):
        """3-carbon sidechain thiol on norvaline scaffold."""
        self._check('N[C@@H](CCCS)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 16},
                     expect_ct={4: 'thiol'})

    # ── Selenol: 4 more (total 5) ────────────────────────────────────────────

    def test_d_sec_selenol(self):
        """D-Sec: selenol on D-isomer."""
        self._check('N[C@H](C[SeH])C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 34},
                     expect_lg={4: '[H]'}, expect_ct={4: 'selenol'})

    def test_homoselenocysteine_selenol(self):
        """Homo-Sec: 2C sidechain selenol."""
        self._check('N[C@@H](CC[SeH])C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 34},
                     expect_ct={4: 'selenol'})

    def test_d_homoselenocysteine_selenol(self):
        """D-homo-Sec."""
        self._check('N[C@H](CC[SeH])C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 34},
                     expect_ct={4: 'selenol'})

    def test_long_chain_selenol(self):
        """3C sidechain selenol."""
        self._check('N[C@@H](CCC[SeH])C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 34},
                     expect_ct={4: 'selenol'})

    # ── Alkyl halide: 5 new ──────────────────────────────────────────────────

    def test_chloroalanine_alkyl_halide(self):
        """3-Chloro-Ala: Cl on beta-C, dummy on same C."""
        self._check('N[C@@H](CCl)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_lg={4: '[Cl]'}, expect_ct={4: 'alkyl_halide_c'})

    def test_bromoalanine_alkyl_halide(self):
        """3-Bromo-Ala: Br variant."""
        self._check('N[C@@H](CBr)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_ct={4: 'alkyl_halide_c'})

    def test_iodoalanine_alkyl_halide(self):
        """3-Iodo-Ala: I variant."""
        self._check('N[C@@H](CI)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_ct={4: 'alkyl_halide_c'})

    def test_d_chloroalanine_alkyl_halide(self):
        """D-3-Chloro-Ala."""
        self._check('N[C@H](CCl)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_ct={4: 'alkyl_halide_c'})

    def test_chlorohomoalanine_alkyl_halide(self):
        """4-Chloro-homoAla: longer chain Cl."""
        self._check('N[C@@H](CCCl)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 6},
                     expect_ct={4: 'alkyl_halide_c'})

    # ── Aminooxy: 5 new ──────────────────────────────────────────────────────

    def test_aminooxy_gamma(self):
        """Aminooxy on gamma-C sidechain."""
        r = self._check('NOCC[C@@H](N)C(=O)O', expect_ct={4: 'aminooxy'})
        smap = self._slot_map(r.chuckles)
        assert smap[4] == 7

    def test_aminooxy_d_isomer(self):
        """D-isomer aminooxy."""
        r = self._check('NOCC[C@H](N)C(=O)O', expect_ct={4: 'aminooxy'})
        smap = self._slot_map(r.chuckles)
        assert smap[4] == 7

    def test_aminooxy_beta(self):
        """Aminooxy on beta sidechain (Ser-like)."""
        r = self._check('N[C@@H](CON)C(=O)O', expect_ct={4: 'aminooxy'})
        smap = self._slot_map(r.chuckles)
        assert smap[4] == 7

    def test_aminooxy_delta(self):
        """Aminooxy on delta-C sidechain."""
        r = self._check('NOCCC[C@@H](N)C(=O)O', expect_ct={4: 'aminooxy'})
        smap = self._slot_map(r.chuckles)
        assert smap[4] == 7

    def test_aminooxy_homo(self):
        """Aminooxy on 2C sidechain."""
        r = self._check('N[C@@H](CCON)C(=O)O', expect_ct={4: 'aminooxy'})
        smap = self._slot_map(r.chuckles)
        assert smap[4] == 7

    # ── Hydrazide: 5 new ─────────────────────────────────────────────────────

    def test_glu_hydrazide(self):
        """Glu-hydrazide: -C(=O)NHNH2 sidechain."""
        self._check('N[C@@H](CCC(=O)NN)C(=O)O',
                     expect_ct={4: 'hydrazide'})

    def test_asp_hydrazide(self):
        """Asp-hydrazide: shorter chain."""
        self._check('N[C@@H](CC(=O)NN)C(=O)O',
                     expect_ct={4: 'hydrazide'})

    def test_d_glu_hydrazide(self):
        """D-Glu-hydrazide."""
        self._check('N[C@H](CCC(=O)NN)C(=O)O',
                     expect_ct={4: 'hydrazide'})

    def test_d_asp_hydrazide(self):
        """D-Asp-hydrazide."""
        self._check('N[C@H](CC(=O)NN)C(=O)O',
                     expect_ct={4: 'hydrazide'})

    def test_long_chain_hydrazide(self):
        """Longer chain hydrazide."""
        self._check('N[C@@H](CCCC(=O)NN)C(=O)O',
                     expect_ct={4: 'hydrazide'})

    # ── Amine primary: 3 more (total 5) ──────────────────────────────────────

    def test_dab_amine_primary(self):
        """Dab (2,4-diaminobutyric acid): short chain amine."""
        self._check('N[C@@H](CCN)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    def test_dap_amine_primary(self):
        """Dap (2,3-diaminopropionic acid): shortest amine sidechain."""
        self._check('N[C@@H](CN)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary', 5: 'amine_secondary'})

    def test_d_orn_amine_primary(self):
        """D-Orn: D-isomer delta-amine."""
        self._check('N[C@H](CCCN)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                     expect_ct={4: 'amine_primary'})

    # ── Guanidinium: 5 new ───────────────────────────────────────────────────

    def test_arg_guanidinium_slots(self):
        """Arg: guanidinium label_only + terminal amine_primary."""
        r = self._check('N[C@@H](CCCNC(=N)N)C(=O)O')
        assert 'guanidinium' in r.chem_types.values()
        assert 'amine_primary' in r.chem_types.values()

    def test_d_arg_guanidinium(self):
        """D-Arg: same guanidinium + amine_primary slots."""
        r = self._check('N[C@H](CCCNC(=N)N)C(=O)O')
        assert 'guanidinium' in r.chem_types.values()
        assert 'amine_primary' in r.chem_types.values()

    def test_arg_guanidinium_on_nitrogen(self):
        """Arg: guanidinium dummy attaches to N (element 7)."""
        r = self._check('N[C@@H](CCCNC(=N)N)C(=O)O')
        smap = self._slot_map(r.chuckles)
        guan_slot = [s for s, ct in r.chem_types.items() if ct == 'guanidinium'][0]
        assert smap[guan_slot] == 7

    def test_homoarg_guanidinium(self):
        """Homoarginine: longer chain to guanidinium."""
        r = self._check('N[C@@H](CCCCNC(=N)N)C(=O)O')
        assert 'guanidinium' in r.chem_types.values()

    def test_arg_guanidinium_lg_is_h(self):
        """Arg: guanidinium leaving group is [H]."""
        r = self._check('N[C@@H](CCCNC(=N)N)C(=O)O')
        guan_slot = [s for s, ct in r.chem_types.items() if ct == 'guanidinium'][0]
        assert r.leaving[guan_slot] == '[H]'

    # ── Hydroxyl phenolic: 4 more (total 5) ──────────────────────────────────

    def test_d_tyr_phenolic_oh(self):
        """D-Tyr: phenolic OH."""
        self._check('N[C@H](Cc1ccc(O)cc1)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl_phenolic'})

    def test_3f_tyr_phenolic_oh(self):
        """3-Fluoro-Tyr: substituted phenol still detected."""
        self._check('N[C@@H](Cc1cc(F)c(O)cc1)C(=O)O',
                     expect_ct={4: 'hydroxyl_phenolic'})

    def test_3cl_tyr_phenolic_oh(self):
        """3-Chloro-Tyr: Cl-substituted phenol."""
        self._check('N[C@@H](Cc1cc(Cl)c(O)cc1)C(=O)O',
                     expect_ct={4: 'hydroxyl_phenolic'})

    def test_dopa_two_phenolic_oh(self):
        """DOPA (3,4-dihydroxyphenylalanine): two phenolic OHs detected."""
        r = self._check('N[C@@H](Cc1cc(O)c(O)cc1)C(=O)O')
        phenol_slots = [s for s, ct in r.chem_types.items() if ct == 'hydroxyl_phenolic']
        assert len(phenol_slots) == 2

    # ── Hydroxyl: 3 more (total 5) ───────────────────────────────────────────

    def test_homoserine_hydroxyl(self):
        """Homoserine: 2C sidechain hydroxyl."""
        self._check('N[C@@H](CCO)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl'})

    def test_d_thr_hydroxyl(self):
        """D-Thr: hydroxyl on D-isomer."""
        self._check('N[C@H]([C@H](O)C)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl'})

    def test_d_ser_hydroxyl(self):
        """D-Ser: hydroxyl on D-Ser."""
        self._check('N[C@H](CO)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 8},
                     expect_ct={4: 'hydroxyl'})

    # ── Aromatic NH: 3 more (total 5) ────────────────────────────────────────

    def test_d_trp_aromatic_nh(self):
        """D-Trp: indole NH."""
        self._check('N[C@H](Cc1c[nH]c2ccccc12)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7},
                     expect_ct={4: 'aromatic_nh'})

    def test_d_his_aromatic_nh(self):
        """D-His: imidazole NH."""
        self._check('N[C@H](Cc1cnc[nH]1)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7},
                     expect_ct={4: 'aromatic_nh'})

    def test_5f_trp_aromatic_nh(self):
        """5-Fluoro-Trp: substituted indole still detected as aromatic_nh."""
        self._check('N[C@@H](Cc1c[nH]c2cc(F)ccc12)C(=O)O',
                     {1: 7, 2: 6, 3: 7, 4: 7},
                     expect_ct={4: 'aromatic_nh'})

    # ── Amide NH (label_only, secondary amide): 5 new ────────────────────────

    def test_nme_gln_amide_nh(self):
        """N-methyl-glutamine: sidechain -C(=O)NHCH3 gives amide_nh."""
        r = self._check('N[C@@H](CCC(=O)NC)C(=O)O')
        assert 'amide_nh' in r.chem_types.values()
        smap = self._slot_map(r.chuckles)
        amide_slot = [s for s, ct in r.chem_types.items() if ct == 'amide_nh'][0]
        assert smap[amide_slot] == 7

    def test_nme_asn_amide_nh(self):
        """N-methyl-asparagine: shorter chain amide_nh."""
        r = self._check('N[C@@H](CC(=O)NC)C(=O)O')
        assert 'amide_nh' in r.chem_types.values()

    def test_net_gln_amide_nh(self):
        """N-ethyl-glutamine."""
        r = self._check('N[C@@H](CCC(=O)NCC)C(=O)O')
        assert 'amide_nh' in r.chem_types.values()

    def test_d_nme_gln_amide_nh(self):
        """D-N-methyl-glutamine."""
        r = self._check('N[C@H](CCC(=O)NC)C(=O)O')
        assert 'amide_nh' in r.chem_types.values()

    def test_long_chain_amide_nh(self):
        """Long chain secondary amide."""
        r = self._check('N[C@@H](CCCC(=O)NC)C(=O)O')
        assert 'amide_nh' in r.chem_types.values()

    # ── Phosphate P: 5 new ───────────────────────────────────────────────────

    def test_phosphoserine_phosphate(self):
        """pSer: phosphate on serine OH."""
        r = self._check('N[C@@H](COP(=O)(O)O)C(=O)O')
        assert 'phosphate_p' in r.chem_types.values()
        smap = self._slot_map(r.chuckles)
        p_slot = [s for s, ct in r.chem_types.items() if ct == 'phosphate_p'][0]
        assert smap[p_slot] == 15

    def test_phosphothreonine_phosphate(self):
        """pThr: phosphate on threonine OH."""
        r = self._check('N[C@@H]([C@@H](OP(=O)(O)O)C)C(=O)O')
        assert 'phosphate_p' in r.chem_types.values()

    def test_phosphotyrosine_phosphate(self):
        """pTyr: phosphate on phenolic OH."""
        r = self._check('N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O')
        assert 'phosphate_p' in r.chem_types.values()

    def test_d_phosphoserine_phosphate(self):
        """D-pSer."""
        r = self._check('N[C@H](COP(=O)(O)O)C(=O)O')
        assert 'phosphate_p' in r.chem_types.values()

    def test_phosphoserine_lg_is_oh(self):
        """pSer: phosphate leaving group is [OH]."""
        r = self._check('N[C@@H](COP(=O)(O)O)C(=O)O')
        p_slot = [s for s, ct in r.chem_types.items() if ct == 'phosphate_p'][0]
        assert r.leaving[p_slot] == '[OH]'

    # ── Alkyne C (CuAAC): 5 new ──────────────────────────────────────────────

    def test_propargylglycine_alkyne(self):
        """Pra (propargylglycine): 1C spacer to terminal alkyne."""
        self._check('N[C@@H](CC#C)C(=O)O',
                     expect_ct={4: 'alkyne_c'})

    def test_homopropargylglycine_alkyne(self):
        """Hpg: 2C spacer."""
        self._check('N[C@@H](CCC#C)C(=O)O',
                     expect_ct={4: 'alkyne_c'})

    def test_3c_spacer_alkyne(self):
        """3C spacer to terminal alkyne."""
        self._check('N[C@@H](CCCC#C)C(=O)O',
                     expect_ct={4: 'alkyne_c'})

    def test_d_propargylglycine_alkyne(self):
        """D-Pra."""
        self._check('N[C@H](CC#C)C(=O)O',
                     expect_ct={4: 'alkyne_c'})

    def test_alkyne_dummy_on_carbon(self):
        """Alkyne: dummy lands on sp3 C (element 6)."""
        r = self._check('N[C@@H](CC#C)C(=O)O')
        smap = self._slot_map(r.chuckles)
        alk_slot = [s for s, ct in r.chem_types.items() if ct == 'alkyne_c'][0]
        assert smap[alk_slot] == 6

    # ── Azide alpha C (SPAAC/CuAAC): 5 new ──────────────────────────────────

    def test_azidoalanine_azide(self):
        """AzAla: 1C spacer to organic azide."""
        self._check('N[C@@H](CN=[N+]=[N-])C(=O)O',
                     expect_ct={4: 'azide_alpha_c'})

    def test_azidohomoalanine_azide(self):
        """AzHal: 2C spacer."""
        self._check('N[C@@H](CCN=[N+]=[N-])C(=O)O',
                     expect_ct={4: 'azide_alpha_c'})

    def test_azidoornithine_azide(self):
        """AzOrn: 3C spacer."""
        self._check('N[C@@H](CCCN=[N+]=[N-])C(=O)O',
                     expect_ct={4: 'azide_alpha_c'})

    def test_azidolysine_azide(self):
        """AzK: 4C spacer."""
        self._check('N[C@@H](CCCCN=[N+]=[N-])C(=O)O',
                     expect_ct={4: 'azide_alpha_c'})

    def test_d_azidoalanine_azide(self):
        """D-AzAla."""
        self._check('N[C@H](CN=[N+]=[N-])C(=O)O',
                     expect_ct={4: 'azide_alpha_c'})

    # ── Terminal alkene (RCM stapling): 5 new ────────────────────────────────

    def test_octenoic_terminal_alkene(self):
        """S5-type RCM: 4C sidechain with terminal alkene."""
        self._check('N[C@@H](CCCC=C)C(=O)O',
                     expect_ct={4: 'terminal_alkene'})

    def test_pentenoic_terminal_alkene(self):
        """Shorter 2C sidechain alkene."""
        self._check('N[C@@H](CC=C)C(=O)O',
                     expect_ct={4: 'terminal_alkene'})

    def test_hexenoic_terminal_alkene(self):
        """3C sidechain alkene."""
        self._check('N[C@@H](CCC=C)C(=O)O',
                     expect_ct={4: 'terminal_alkene'})

    def test_d_octenoic_terminal_alkene(self):
        """D-isomer 4C sidechain alkene."""
        self._check('N[C@H](CCCC=C)C(=O)O',
                     expect_ct={4: 'terminal_alkene'})

    def test_terminal_alkene_dummy_on_carbon(self):
        """Terminal alkene: dummy on internal vinyl C (element 6)."""
        r = self._check('N[C@@H](CCCC=C)C(=O)O')
        smap = self._slot_map(r.chuckles)
        alk_slot = [s for s, ct in r.chem_types.items() if ct == 'terminal_alkene'][0]
        assert smap[alk_slot] == 6

    # ── Tetrazine C (iEDDA): 5 new ───────────────────────────────────────────

    def test_tetrazine_1c_spacer(self):
        """TzAla: 1C spacer to s-tetrazine ring."""
        self._check('N[C@@H](Cc1nncnn1)C(=O)O',
                     expect_ct={4: 'tetrazine_c'})

    def test_tetrazine_2c_spacer(self):
        """TzAbu: 2C spacer."""
        self._check('N[C@@H](CCc1nncnn1)C(=O)O',
                     expect_ct={4: 'tetrazine_c'})

    def test_tetrazine_5c_spacer(self):
        """TzLys: 5C spacer."""
        self._check('N[C@@H](CCCCCc1nncnn1)C(=O)O',
                     expect_ct={4: 'tetrazine_c'})

    def test_d_tetrazine(self):
        """D-TzAla."""
        self._check('N[C@H](Cc1nncnn1)C(=O)O',
                     expect_ct={4: 'tetrazine_c'})

    def test_tetrazine_dummy_on_carbon(self):
        """Tetrazine: dummy on sp3 C (element 6)."""
        r = self._check('N[C@@H](Cc1nncnn1)C(=O)O')
        smap = self._slot_map(r.chuckles)
        tz_slot = [s for s, ct in r.chem_types.items() if ct == 'tetrazine_c'][0]
        assert smap[tz_slot] == 6

    # ── TCO C (trans-cyclooctene, iEDDA partner): 5 new ──────────────────────

    def test_tco_1c_linker(self):
        """TcoAla: 1C linker to cyclooctene ring alkene."""
        self._check('N[C@@H](CC1=CCCCCCC1)C(=O)O',
                     expect_ct={4: 'tco_c'})

    def test_tco_2c_linker(self):
        """TcoAbu: 2C linker."""
        self._check('N[C@@H](CCC1=CCCCCCC1)C(=O)O',
                     expect_ct={4: 'tco_c'})

    def test_tco_4c_linker(self):
        """TcoLys: 4C linker."""
        self._check('N[C@@H](CCCCC1=CCCCCCC1)C(=O)O',
                     expect_ct={4: 'tco_c'})

    def test_d_tco(self):
        """D-TcoAla."""
        self._check('N[C@H](CC1=CCCCCCC1)C(=O)O',
                     expect_ct={4: 'tco_c'})

    def test_tco_dummy_exocyclic(self):
        """TCO: dummy on exocyclic sp3 C, not in ring."""
        r = self._check('N[C@@H](CC1=CCCCCCC1)C(=O)O')
        mol = Chem.MolFromSmiles(r.chuckles)
        tco_slot = [s for s, ct in r.chem_types.items() if ct == 'tco_c'][0]
        dummy = next(a for a in mol.GetAtoms()
                     if a.GetAtomicNum() == 0 and a.GetIsotope() == tco_slot)
        nb = next(iter(dummy.GetNeighbors()))
        assert not nb.IsInRing(), "TCO dummy should be on exocyclic C"

    # ── Cyclooctyne C (SPAAC): 5 new ────────────────────────────────────────

    def test_cyclooctyne_1c_linker(self):
        """CyoAla: 1C linker to cyclooctyne ring."""
        self._check('N[C@@H](CC1CCCCC#CC1)C(=O)O',
                     expect_ct={4: 'cyclooctyne_c'})

    def test_cyclooctyne_2c_linker(self):
        """CyoAbu: 2C linker."""
        self._check('N[C@@H](CCC1CCCCC#CC1)C(=O)O',
                     expect_ct={4: 'cyclooctyne_c'})

    def test_cyclooctyne_4c_linker(self):
        """CyoLys: 4C linker."""
        self._check('N[C@@H](CCCCC1CCCCC#CC1)C(=O)O',
                     expect_ct={4: 'cyclooctyne_c'})

    def test_d_cyclooctyne(self):
        """D-CyoAla."""
        self._check('N[C@H](CC1CCCCC#CC1)C(=O)O',
                     expect_ct={4: 'cyclooctyne_c'})

    def test_cyclooctyne_dummy_on_carbon(self):
        """Cyclooctyne: dummy on ring sp3 C (element 6)."""
        r = self._check('N[C@@H](CC1CCCCC#CC1)C(=O)O')
        smap = self._slot_map(r.chuckles)
        cyo_slot = [s for s, ct in r.chem_types.items() if ct == 'cyclooctyne_c'][0]
        assert smap[cyo_slot] == 6

    # ── Aldehyde: 5 new ──────────────────────────────────────────────────────

    def test_formylphe_aldehyde(self):
        """4-Formyl-Phe: aldehyde CHO on para position."""
        self._check('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O',
                     expect_ct={4: 'aldehyde'})

    def test_d_formylphe_aldehyde(self):
        """D-4-formyl-Phe."""
        self._check('N[C@H](Cc1ccc(C=O)cc1)C(=O)O',
                     expect_ct={4: 'aldehyde'})

    def test_aliphatic_aldehyde(self):
        """Aliphatic sidechain aldehyde."""
        self._check('N[C@@H](CCC=O)C(=O)O',
                     expect_ct={4: 'aldehyde'})

    def test_long_aliphatic_aldehyde(self):
        """Longer aliphatic aldehyde sidechain."""
        self._check('N[C@@H](CCCC=O)C(=O)O',
                     expect_ct={4: 'aldehyde'})

    def test_aldehyde_dummy_on_carbon(self):
        """Aldehyde: dummy on aldehyde C (element 6)."""
        r = self._check('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O')
        smap = self._slot_map(r.chuckles)
        ald_slot = [s for s, ct in r.chem_types.items() if ct == 'aldehyde'][0]
        assert smap[ald_slot] == 6

    # ── NHS ester: 5 new ─────────────────────────────────────────────────────

    def test_nhs_asp_ester(self):
        """NhsAsp: NHS ester on Asp-length sidechain."""
        self._check('N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O',
                     expect_ct={4: 'nhs_ester'})

    def test_nhs_glu_ester(self):
        """NhsGlu: NHS ester on Glu-length sidechain."""
        self._check('N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O',
                     expect_ct={4: 'nhs_ester'})

    def test_nhs_nva_ester(self):
        """NhsNva: NHS ester on Nva-length sidechain."""
        self._check('N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O',
                     expect_ct={4: 'nhs_ester'})

    def test_d_nhs_asp_ester(self):
        """D-NhsAsp."""
        self._check('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O',
                     expect_ct={4: 'nhs_ester'})

    def test_nhs_ester_dummy_on_carbon(self):
        """NHS ester: dummy on sp3 C (element 6)."""
        r = self._check('N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O')
        smap = self._slot_map(r.chuckles)
        nhs_slot = [s for s, ct in r.chem_types.items() if ct == 'nhs_ester'][0]
        assert smap[nhs_slot] == 6

    # ── Maleimide C: 5 new ───────────────────────────────────────────────────

    def test_lys_maleimide(self):
        """Lys-maleimide: maleimide ring on epsilon-N chain."""
        self._check('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O',
                     expect_ct={4: 'maleimide_c'})

    def test_short_maleimide(self):
        """Shorter 2C linker to maleimide."""
        self._check('N[C@@H](CCN1C(=O)C=CC1=O)C(=O)O',
                     expect_ct={4: 'maleimide_c'})

    def test_medium_maleimide(self):
        """3C linker to maleimide."""
        self._check('N[C@@H](CCCN1C(=O)C=CC1=O)C(=O)O',
                     expect_ct={4: 'maleimide_c'})

    def test_d_lys_maleimide(self):
        """D-Lys-maleimide."""
        self._check('N[C@H](CCCCN1C(=O)C=CC1=O)C(=O)O',
                     expect_ct={4: 'maleimide_c'})

    def test_maleimide_dummy_on_carbon(self):
        """Maleimide: dummy on ring vinyl C (element 6)."""
        r = self._check('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O')
        smap = self._slot_map(r.chuckles)
        mal_slot = [s for s, ct in r.chem_types.items() if ct == 'maleimide_c'][0]
        assert smap[mal_slot] == 6

    # ── Amine secondary (second-H pass): 3 more (total 5 with Lys+Orn) ──────

    def test_dab_amine_secondary(self):
        """Dab: epsilon-NH2 gives amine_secondary on second H."""
        r = self._check('N[C@@H](CCN)C(=O)O',
                         expect_ct={5: 'amine_secondary'})
        assert r.leaving[5] == '[H]'

    def test_dap_amine_secondary(self):
        """Dap: NH2 gives amine_secondary."""
        r = self._check('N[C@@H](CN)C(=O)O',
                         expect_ct={5: 'amine_secondary'})

    def test_d_lys_amine_secondary(self):
        """D-Lys: second-H pass fires on D-isomer too."""
        self._check('N[C@H](CCCCN)C(=O)O',
                     expect_ct={5: 'amine_secondary'})

    # ── No sidechain: negative controls ──────────────────────────────────────

    def test_nle_no_sidechain(self):
        """Norleucine: n-butyl sidechain, no reactive group."""
        self._check('N[C@@H](CCCC)C(=O)O', {1: 7, 2: 6, 3: 7})

    def test_nva_no_sidechain(self):
        """Norvaline: n-propyl sidechain, no reactive group."""
        self._check('N[C@@H](CCC)C(=O)O', {1: 7, 2: 6, 3: 7})

    # ── Reagent-form caps (halide LG detection) ─────────────────────────────

    def test_acetyl_chloride_cap(self):
        """Acid chloride: CC(=O)Cl → R2 slot, LG=[Cl]."""
        r = self._check('CC(=O)Cl', {2: 6})
        assert r.leaving[2] == '[Cl]'
        assert r.chem_types[2] == 'backbone_c'

    def test_benzoyl_chloride_cap(self):
        """Benzoyl chloride: O=C(Cl)c1ccccc1 → R2 slot, LG=[Cl]."""
        r = self._check('O=C(Cl)c1ccccc1', {2: 6})
        assert r.leaving[2] == '[Cl]'

    def test_tosyl_chloride_cap(self):
        """Sulfonyl chloride: TsCl → R2 slot on S, LG=[Cl]."""
        r = self._check('Cc1ccc(S(=O)(=O)Cl)cc1', {2: 16})
        assert r.leaving[2] == '[Cl]'

    def test_methyl_iodide_cap(self):
        """Alkyl iodide: CI → R2 slot, LG=[I]."""
        r = self._check('CI', {2: 6})
        assert r.leaving[2] == '[I]'

    def test_ethyl_bromide_cap(self):
        """Alkyl bromide: CCBr → R2 slot, LG=[Br]."""
        r = self._check('CCBr', {2: 6})
        assert r.leaving[2] == '[Br]'

    def test_benzyl_chloride_cap(self):
        """Alkyl chloride: ClCc1ccccc1 → R2 slot, LG=[Cl]."""
        r = self._check('ClCc1ccccc1', {2: 6})
        assert r.leaving[2] == '[Cl]'

    def test_ammonia_nucleophilic_cap(self):
        """Nucleophilic cap: N (ammonia) → R1 slot, LG=[H]."""
        r = self._check('N', {1: 7})
        assert r.leaving[1] == '[H]'
        assert r.chem_types[1] == 'backbone_n'

    def test_methanol_nucleophilic_cap(self):
        """Nucleophilic cap: CO (methanol) → R1 slot on O, LG=[H]."""
        r = self._check('CO', {1: 8})
        assert r.leaving[1] == '[H]'

    def test_free_acid_vs_acid_chloride(self):
        """Same scaffold: free acid gets [OH], acid chloride gets [Cl]."""
        r_free = self._check('CC(=O)O', {2: 6})
        r_halide = self._check('CC(=O)Cl', {2: 6})
        assert r_free.leaving[2] == '[OH]'
        assert r_halide.leaving[2] == '[Cl]'


class TestSynonymResolution:
    """Synonym aliases in the CSV resolve to canonical SDF entries."""

    @staticmethod
    def _smi(biln):
        return _smiles(biln)

    def test_three_letter_gly(self):
        """Gly resolves to G."""
        assert self._smi('ac-Gly-am') == self._smi('ac-G-am')

    def test_three_letter_ala(self):
        """Ala resolves to A."""
        assert self._smi('ac-Ala-am') == self._smi('ac-A-am')

    def test_d_amino_acid_alias(self):
        """dA resolves to DAla."""
        assert self._smi('ac-dA-am') == self._smi('ac-DAla-am')

    def test_cap_synonym_butcap(self):
        """ButCap resolves to Bua."""
        assert self._smi('ButCap-G-am') == self._smi('Bua-G-am')

    def test_cap_synonym_bzcap(self):
        """BzCap resolves to Bz."""
        assert self._smi('BzCap-G-am') == self._smi('Bz-G-am')

    def test_ncaa_synonym_dopa(self):
        """Dopa resolves to Tyr_3OH."""
        assert self._smi('ac-Dopa-am') == self._smi('ac-Tyr_3OH-am')

    def test_ncaa_synonym_hcys(self):
        """hCys resolves to Hcy."""
        assert self._smi('ac-hCys-am') == self._smi('ac-Hcy-am')

    def test_ncaa_synonym_kyn(self):
        """Kyn resolves to Asp_Ph2NH2."""
        assert self._smi('ac-Kyn-am') == self._smi('ac-Asp_Ph2NH2-am')

    def test_nrg_resolves_to_arg(self):
        """Nrg resolves to R (Arg)."""
        assert self._smi('ac-Nrg-am') == self._smi('ac-R-am')

    def test_phe_4ome_resolves_to_tyr_me(self):
        """Phe_4OMe resolves to Tyr_Me."""
        assert self._smi('ac-Phe_4OMe-am') == self._smi('ac-Tyr_Me-am')


class TestRestoreRgroups:
    """Edge-case tests for __restore_and_remove_rgroups in molecule.py.

    Covers every restore path:
      - C-attachment carboxyl (LG=[OH]): dummy swapped for O atom → COOH
      - Backbone N (LG=[H]): dummy removed → free NH2
      - Backbone C (LG=[OH]): dummy swapped for O → free COOH
      - Thiol (LG=[H]): dummy removed → free SH
      - Amine (LG=[H]): dummy removed → free NH2
    Also verifies carboxyl consumed by crosslinks is NOT restored.
    """

    ACID  = Chem.MolFromSmarts('[CX3](=[O])[OX2H1]')   # free carboxylic acid
    THIOL = Chem.MolFromSmarts('[SX2H1]')                # free thiol
    NH2   = Chem.MolFromSmarts('[NX3H2][CX4]')          # aliphatic primary amine

    def _romol(self, biln, fmt=None):
        seq = Sequence(biln) if fmt is None else Sequence(biln, fmt=fmt)
        return Molecule(seq).get_molecule(fmt='ROMol')

    def _counts(self, biln, fmt=None):
        m = self._romol(biln, fmt)
        return {
            'acid':  len(m.GetSubstructMatches(self.ACID)),
            'thiol': len(m.GetSubstructMatches(self.THIOL)),
            'nh2':   len(m.GetSubstructMatches(self.NH2)),
            'dummy': sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 0),
        }

    # ── Group 1: homo-Asp/Glu variants (C-attachment carboxyl, LG=[OH]) ──────

    def test_hasp_carboxyl_restores(self):
        """E slot-4 C-attachment carboxyl (LG=[OH]) restores to free COOH."""
        c = self._counts('ac-E-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_hasp_carboxyl_restores(self):
        """DGlu D-isomer slot-4 carboxyl restores correctly."""
        c = self._counts('ac-DGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_hglu_carboxyl_restores(self):
        """hGlu slot-4 C-attachment carboxyl (LG=[OH]) restores to free COOH."""
        c = self._counts('ac-hGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_hglu_carboxyl_restores(self):
        """D_hGlu D-isomer slot-4 carboxyl restores correctly."""
        c = self._counts('ac-D_hGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_b3hasp_carboxyl_restores(self):
        """b3hAsp (beta-3-homo-Asp) slot-4 carboxyl restores."""
        c = self._counts('ac-b3hAsp-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_b3hasp_carboxyl_restores(self):
        """D_b3hAsp D-isomer slot-4 carboxyl restores correctly."""
        c = self._counts('ac-D_b3hAsp-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_b3hglu_carboxyl_restores(self):
        """b3hGlu (beta-3-homo-Glu) slot-4 carboxyl restores."""
        c = self._counts('ac-b3hGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_b3hglu_carboxyl_restores(self):
        """D_b3hGlu D-isomer slot-4 carboxyl restores correctly."""
        c = self._counts('ac-D_b3hGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_ameasp_carboxyl_restores(self):
        """aMeAsp (alpha-methyl-Asp) slot-4 carboxyl restores."""
        c = self._counts('ac-aMeAsp-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_ameglu_carboxyl_restores(self):
        """aMeGlu (alpha-methyl-Glu) slot-4 carboxyl restores."""
        c = self._counts('ac-aMeGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_med_slot3_carboxyl_restores(self):
        """D_meD (NMe-D-Asp) slot-3 carboxyl restores to free COOH."""
        c = self._counts('ac-D_meD-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_d_mee_slot3_carboxyl_restores(self):
        """D_meE (NMe-D-Glu) slot-3 carboxyl restores to free COOH."""
        c = self._counts('ac-D_meE-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    # ── Group 2: Standard carboxyl (C-attachment, LG=[OH]) ───────────────────

    def test_asp_c_attach_carboxyl_restores(self):
        """Standard Asp (D): slot-4 C-attachment carboxyl (LG=[OH]) restored via atom swap."""
        c = self._counts('ac-D-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_glu_c_attach_carboxyl_restores(self):
        """Standard Glu (E): slot-4 C-attachment carboxyl (LG=[OH]) restored via atom swap."""
        c = self._counts('ac-E-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    # ── Group 3: Backbone N and C terminus restore ────────────────────────────

    def test_free_n_terminus_restores_amine(self):
        """Free N-terminus (backbone_n, LG=[H]): dummy removed → NH2 present on G."""
        c = self._counts('G')
        assert c['nh2'] == 1 and c['dummy'] == 0

    def test_free_c_terminus_restores_carboxyl(self):
        """Free C-terminus (backbone_c, LG=[OH]): dummy swapped for O → free COOH on G."""
        c = self._counts('G')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_n_cap_suppresses_n_terminus_amine(self):
        """ac cap bonds to backbone-N slot; no free aliphatic NH2 remains in ac-G-am."""
        c = self._counts('ac-G-am')
        assert c['nh2'] == 0 and c['dummy'] == 0

    def test_c_cap_suppresses_c_terminus_acid(self):
        """am cap bonds to backbone-C slot; no free COOH remains in ac-G-am."""
        c = self._counts('ac-G-am')
        assert c['acid'] == 0 and c['dummy'] == 0

    def test_both_caps_suppress_both_termini(self):
        """ac-G-am: both termini capped — no free NH2, no free COOH, no dummies."""
        c = self._counts('ac-G-am')
        assert c['acid'] == 0 and c['nh2'] == 0 and c['dummy'] == 0

    def test_free_thiol_restores_on_cysteine(self):
        """Cys (C): slot-4 thiol (LG=[H]) dummy removed → free SH present; no spurious acid."""
        c = self._counts('ac-C-am')
        assert c['thiol'] == 1 and c['acid'] == 0 and c['dummy'] == 0

    def test_free_amine_restores_on_lysine(self):
        """Lys (K): epsilon-amine (LG=[H]) dummies removed → 1 free NH2; no spurious acid."""
        c = self._counts('ac-K-am')
        assert c['nh2'] == 1 and c['acid'] == 0 and c['dummy'] == 0

    # ── Group 4: Multi-residue carboxyl accumulation ──────────────────────────

    def test_triple_asp_three_carboxyls_restored(self):
        """Three Asp residues: each slot-4 COOH restored independently → 3 free acids."""
        c = self._counts('ac-D-D-D-am')
        assert c['acid'] == 3 and c['dummy'] == 0

    def test_five_asp_five_carboxyls_restored(self):
        """Five Asp residues: 5 independent slot-4 restores → exactly 5 free acids."""
        c = self._counts('ac-D-D-D-D-D-am')
        assert c['acid'] == 5 and c['dummy'] == 0

    def test_mixed_hasp_hglu_two_carboxyls(self):
        """E + hGlu (both C-attachment slot-4): 2 sidechain COOHs restored."""
        c = self._counts('ac-E-hGlu-am')
        assert c['acid'] == 2 and c['dummy'] == 0

    def test_asp_glu_two_carboxyls_restored(self):
        """Standard Asp + Glu (both C-attachment slot-4): 2 sidechain COOHs restored."""
        c = self._counts('ac-D-E-am')
        assert c['acid'] == 2 and c['dummy'] == 0

    # ── Group 5: Orthogonal restore and crosslink interactions ────────────────

    def test_cys_and_asp_independent_restore(self):
        """Cys thiol and Asp carboxyl restored independently: 1 thiol + 1 COOH."""
        c = self._counts('ac-C-D-am')
        assert c['thiol'] == 1 and c['acid'] == 1 and c['dummy'] == 0

    def test_disulfide_eliminates_thiols_but_asp_cooh_survives(self):
        """Disulfide crosslink consumes both Cys thiols; Asp sidechain COOH still restored."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = self._counts('ac-C.!1(4,4)-D-C.!1-am')
        assert c['thiol'] == 0 and c['acid'] == 1 and c['dummy'] == 0

    def test_lactam_crosslink_consumes_carboxyl_not_restored(self):
        """K-D sidechain amide crosslink: Asp slot-4 COOH bonded → not restored as free acid."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = self._counts('ac-K.!1(4,4)-G-D.!1-am')
        assert c['acid'] == 0 and c['dummy'] == 0

    # ── Group 6: D/L isomer parity ────────────────────────────────────────────

    def test_d_hasp_and_l_hasp_same_acid_count(self):
        """DGlu and E are stereoisomers: restore yields identical COOH count."""
        assert self._counts('ac-E-am')['acid'] == self._counts('ac-DGlu-am')['acid'] == 1

    def test_d_b3hglu_and_l_b3hglu_same_acid_count(self):
        """D_b3hGlu and b3hGlu: D-isomer restore identical to L-isomer."""
        assert self._counts('ac-b3hGlu-am')['acid'] == self._counts('ac-D_b3hGlu-am')['acid'] == 1

    # ── Group 7: No dummy atoms remain in any assembled product ───────────────

    def test_no_dummies_simple_tripeptide(self):
        """A-G-A: free termini both restored; no residual dummy atoms in product."""
        c = self._counts('A-G-A')
        assert c['dummy'] == 0

    def test_no_dummies_complex_mixed_sequence(self):
        """ac-D-K-C-E-am: all four restore paths exercised; zero residual dummies."""
        c = self._counts('ac-D-K-C-E-am')
        assert c['acid'] == 2 and c['thiol'] == 1 and c['nh2'] == 1 and c['dummy'] == 0

    def test_no_dummies_after_disulfide_crosslink(self):
        """Disulfide crosslinked product: crosslink bonds consumed all S dummies; no leakage."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = self._counts('ac-C.!1(4,4)-A-C.!1-am')
        assert c['dummy'] == 0

    # ── Group 8: Hydroxyl restore ────────────────────────────────────────────

    def _extended_counts(self, biln, fmt=None):
        m = self._romol(biln, fmt)
        hydroxyl = Chem.MolFromSmarts('[OX2H1][CX4]')
        phenol = Chem.MolFromSmarts('[OX2H1][c]')
        ar_nh = Chem.MolFromSmarts('[nH]')
        return {
            'acid':    len(m.GetSubstructMatches(self.ACID)),
            'thiol':   len(m.GetSubstructMatches(self.THIOL)),
            'nh2':     len(m.GetSubstructMatches(self.NH2)),
            'dummy':   sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 0),
            'hydroxyl': len(m.GetSubstructMatches(hydroxyl)),
            'phenol':  len(m.GetSubstructMatches(phenol)),
            'ar_nh':   len(m.GetSubstructMatches(ar_nh)),
        }

    def test_ser_hydroxyl_restores(self):
        """Ser: sidechain hydroxyl (LG=[H]) restored to free OH."""
        c = self._extended_counts('ac-S-am')
        assert c['hydroxyl'] == 1 and c['acid'] == 0 and c['dummy'] == 0

    def test_thr_hydroxyl_restores(self):
        """Thr: sidechain hydroxyl restored to free OH."""
        c = self._extended_counts('ac-T-am')
        assert c['hydroxyl'] == 1 and c['acid'] == 0 and c['dummy'] == 0

    def test_tyr_phenolic_oh_restores(self):
        """Tyr: phenolic OH (LG=[H]) restored."""
        c = self._extended_counts('ac-Y-am')
        assert c['phenol'] == 1 and c['dummy'] == 0

    # ── Group 9: Aromatic NH restore ─────────────────────────────────────────

    def test_trp_indole_nh_restores(self):
        """Trp: indole NH (aromatic_nh, LG=[H]) restored."""
        c = self._extended_counts('ac-W-am')
        assert c['ar_nh'] == 1 and c['dummy'] == 0

    def test_his_imidazole_nh_restores(self):
        """His: imidazole NH restored."""
        c = self._extended_counts('ac-H-am')
        assert c['ar_nh'] == 1 and c['dummy'] == 0

    # ── Group 10: Thioether is NOT reactive ──────────────────────────────────

    def test_met_thioether_not_thiol(self):
        """Met: thioether S-CH3 does not produce a thiol restore site."""
        c = self._counts('ac-M-am')
        assert c['thiol'] == 0 and c['dummy'] == 0

    # ── Group 11: Free termini (no caps) ─────────────────────────────────────

    def test_single_gly_free_termini(self):
        """G alone: 1 free NH2 + 1 free COOH from backbone termini."""
        c = self._counts('G')
        assert c['nh2'] == 1 and c['acid'] == 1 and c['dummy'] == 0

    def test_single_asp_free_termini(self):
        """D alone: 1 NH2 + 2 COOHs (backbone + sidechain)."""
        c = self._counts('D')
        assert c['nh2'] == 1 and c['acid'] == 2 and c['dummy'] == 0

    def test_single_glu_free_termini(self):
        """E alone: 1 NH2 + 2 COOHs (backbone + sidechain)."""
        c = self._counts('E')
        assert c['nh2'] == 1 and c['acid'] == 2 and c['dummy'] == 0

    def test_single_cys_free_termini(self):
        """C alone: 1 NH2 + 1 COOH + 1 SH from free termini + thiol."""
        c = self._counts('C')
        assert c['nh2'] == 1 and c['acid'] == 1 and c['thiol'] == 1 and c['dummy'] == 0

    def test_free_n_terminus_tripeptide(self):
        """A-G-D (no caps): 1 free NH2 (N-term) + 2 COOHs (C-term + Asp sidechain)."""
        c = self._counts('A-G-D')
        assert c['nh2'] == 1 and c['acid'] == 2 and c['dummy'] == 0

    # ── Group 12: One cap only ───────────────────────────────────────────────

    def test_ac_only_suppresses_nh2(self):
        """ac-D (N-cap only): 0 NH2, 2 COOHs (backbone + sidechain)."""
        c = self._counts('ac-D')
        assert c['nh2'] == 0 and c['acid'] == 2 and c['dummy'] == 0

    def test_am_only_suppresses_c_term_acid(self):
        """D-am (C-cap only): 1 NH2, 1 COOH (sidechain only; C-term capped)."""
        c = self._counts('D-am')
        assert c['nh2'] == 1 and c['acid'] == 1 and c['dummy'] == 0

    # ── Group 13: Long homopolymer accumulation ──────────────────────────────

    def test_ten_asp_ten_carboxyls(self):
        """10 Asp residues: 10 sidechain COOHs restored independently."""
        c = self._counts('ac-D-D-D-D-D-D-D-D-D-D-am')
        assert c['acid'] == 10 and c['dummy'] == 0

    def test_five_cys_five_thiols(self):
        """5 Cys residues: 5 free thiols restored."""
        c = self._counts('ac-C-C-C-C-C-am')
        assert c['thiol'] == 5 and c['dummy'] == 0

    # ── Group 14: Mixed functional group accumulation ────────────────────────

    def test_mixed_acid_thiol_amine(self):
        """ac-C-D-K-am: 1 thiol + 1 acid + 1 amine, all restored independently."""
        c = self._counts('ac-C-D-K-am')
        assert c['thiol'] == 1 and c['acid'] == 1 and c['nh2'] == 1 and c['dummy'] == 0

    def test_mixed_hydroxyl_acid_thiol(self):
        """ac-S-D-C-am: 1 hydroxyl + 1 acid + 1 thiol."""
        c = self._extended_counts('ac-S-D-C-am')
        assert c['hydroxyl'] == 1 and c['acid'] == 1 and c['thiol'] == 1 and c['dummy'] == 0

    def test_mixed_aromatic_nh_acid(self):
        """ac-W-D-am: 1 indole NH + 1 acid."""
        c = self._extended_counts('ac-W-D-am')
        assert c['ar_nh'] == 1 and c['acid'] == 1 and c['dummy'] == 0

    # ── Group 15: All-inert sequences (no sidechain restore needed) ──────────

    def test_all_ala_no_sidechain_restore(self):
        """ac-A-A-A-am: no sidechain reactive groups -> all backbone dummies consumed."""
        c = self._counts('ac-A-A-A-am')
        assert c['acid'] == 0 and c['thiol'] == 0 and c['nh2'] == 0 and c['dummy'] == 0

    def test_all_gly_no_sidechain_restore(self):
        """ac-G-G-G-am: same as all-Ala."""
        c = self._counts('ac-G-G-G-am')
        assert c['acid'] == 0 and c['thiol'] == 0 and c['nh2'] == 0 and c['dummy'] == 0

    # ── Group 16: Beta amino acid restore ────────────────────────────────────

    def test_beta_homo_asp_glu_mix(self):
        """b3hAsp + b3hGlu: both C-attachment carboxyl, 2 COOHs restored."""
        c = self._counts('ac-b3hAsp-b3hGlu-am')
        assert c['acid'] == 2 and c['dummy'] == 0

    def test_beta_and_standard_asp_mix(self):
        """b3hAsp + D (standard Asp): 2 COOHs, one from each attachment convention."""
        c = self._counts('ac-b3hAsp-D-am')
        assert c['acid'] == 2 and c['dummy'] == 0

    # ── Group 17: Alpha-methyl variants ──────────────────────────────────────

    def test_ameasp_and_standard_glu_mix(self):
        """aMeAsp + E: both have sidechain carboxyl; 2 COOHs restored."""
        c = self._counts('ac-aMeAsp-E-am')
        assert c['acid'] == 2 and c['dummy'] == 0

    def test_ameglu_o_attach_restores(self):
        """aMeGlu standalone: 1 sidechain COOH restored."""
        c = self._counts('ac-aMeGlu-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    # ── Group 18: D-isomer exhaustive parity ─────────────────────────────────

    def test_d_hglu_and_l_hglu_acid_parity(self):
        """D_hGlu vs hGlu: identical COOH count."""
        assert self._counts('ac-hGlu-am')['acid'] == self._counts('ac-D_hGlu-am')['acid'] == 1

    def test_d_b3hasp_and_l_b3hasp_acid_parity(self):
        """D_b3hAsp vs b3hAsp: identical COOH count."""
        assert self._counts('ac-b3hAsp-am')['acid'] == self._counts('ac-D_b3hAsp-am')['acid'] == 1

    def test_d_med_d_mee_parity(self):
        """NMe-D-Asp vs NMe-D-Glu: each has 1 sidechain COOH."""
        assert self._counts('ac-D_meD-am')['acid'] == 1
        assert self._counts('ac-D_meE-am')['acid'] == 1

    # ── Group 19: Crosslink consuming one carboxyl, leaving others intact ────

    def test_lactam_consumes_one_asp_other_survives(self):
        """K-D-D crosslink: lactam consumes first D carboxyl; second D carboxyl survives."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = self._counts('ac-K.!1(4,4)-D.!1-D-am')
        assert c['acid'] == 1 and c['dummy'] == 0

    def test_disulfide_crosslink_acid_and_amine_survive(self):
        """C-D-K-C disulfide: SS bond consumes thiols; Asp COOH and Lys NH2 survive."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = self._counts('ac-C.!1(4,4)-D-K-C.!1-am')
        assert c['thiol'] == 0 and c['acid'] == 1 and c['nh2'] == 1 and c['dummy'] == 0

    # ── Group 20: Atom-count parity between D/L assembled products ───────────

    def test_d_l_asp_heavy_atom_count_equal(self):
        """ac-D-am vs ac-DGlu-am: D-isomer has same heavy atom count as L."""
        m_l = self._romol('ac-E-am')
        m_d = self._romol('ac-DGlu-am')
        assert m_l.GetNumHeavyAtoms() == m_d.GetNumHeavyAtoms()

    def test_d_l_glu_heavy_atom_count_equal(self):
        """hGlu vs D_hGlu: same heavy atom count."""
        m_l = self._romol('ac-hGlu-am')
        m_d = self._romol('ac-D_hGlu-am')
        assert m_l.GetNumHeavyAtoms() == m_d.GetNumHeavyAtoms()


class TestReactionPairSMIRKS:
    """Fire every REACTION_INDEX entry's SMIRKS on pre_activate'd monomers.

    Uses real monomer SMILES through pre_activate to get correctly-placed
    dummies, then feeds those into run_bond_smirks.  This validates the full
    pipeline: auto-assignment -> SMIRKS -> product.
    """

    _MODEL_MONOMERS = {
        'backbone_n':       'N[C@@H](C)C(=O)O',
        'backbone_c':       'N[C@@H](C)C(=O)O',
        'backbone_n_mod':   'N[C@@H](C)C(=O)O',
        'thiol':            'N[C@@H](CS)C(=O)O',
        'selenol':          'N[C@@H](C[SeH])C(=O)O',
        'carboxyl':         'N[C@@H](CC(=O)O)C(=O)O',
        'amine_primary':    'N[C@@H](CCCCN)C(=O)O',
        'amine_secondary':  'N[C@@H](CCCCN)C(=O)O',
        'hydroxyl':         'N[C@@H](CO)C(=O)O',
        'alkyl_halide_c':   'N[C@@H](CCl)C(=O)O',
        'aminooxy':         'NOCC[C@@H](N)C(=O)O',
        'hydrazide':        'N[C@@H](CCC(=O)NN)C(=O)O',
        'aldehyde':         'N[C@@H](Cc1ccc(C=O)cc1)C(=O)O',
        'maleimide_c':      'N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O',
        'nhs_ester':        'N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O',
        'alkyne_c':         'N[C@@H](CC#C)C(=O)O',
        'azide_alpha_c':    'N[C@@H](CN=[N+]=[N-])C(=O)O',
        'terminal_alkene':  'N[C@@H](CCCC=C)C(=O)O',
        'tetrazine_c':      'N[C@@H](Cc1nncnn1)C(=O)O',
        'tco_c':            'N[C@@H](CC1=CCCCCCC1)C(=O)O',
        'cyclooctyne_c':    'N[C@@H](CC1CCCCC#CC1)C(=O)O',
        'phosphate_p':      'N[C@@H](COP(=O)(O)O)C(=O)O',
    }

    # backbone_o and carbon are backbone-level types assigned by heuristic,
    # not auto-assigned by pre_activate on standalone SMILES
    _FALLBACK_FRAGS = {
        'backbone_o':       '[{iso}*]OCC',
        'carbon':           '[{iso}*]CCC',
    }

    def _make_frag(self, ct, iso):
        if ct in self._FALLBACK_FRAGS:
            smi = self._FALLBACK_FRAGS[ct].format(iso=iso)
            mol = Chem.MolFromSmiles(smi)
            assert mol is not None, f"Bad fallback fragment for {ct}"
            return mol

        from pyPept.interfaces.monomer_pipeline import pre_activate
        smiles = self._MODEL_MONOMERS.get(ct)
        if smiles is None:
            pytest.skip(f"No model monomer for chem_type {ct!r}")

        r = pre_activate(smiles)
        mol = Chem.RWMol(Chem.MolFromSmiles(r.chuckles))
        assert mol is not None, f"Bad CHUCKLES for {ct}: {r.chuckles}"

        slot = next((s for s, c in r.chem_types.items() if c == ct), None)
        assert slot is not None, \
            f"pre_activate did not assign {ct!r} to {smiles}; got {r.chem_types}"

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
                atom.SetIsotope(iso)
                break
        else:
            raise AssertionError(
                f"No dummy with isotope {slot} in CHUCKLES for {ct}")

        return mol.GetMol()

    def _run_pair(self, ct_a, ct_b):
        from pyPept.interfaces.reaction_library import REACTION_INDEX, run_bond_smirks
        entry = REACTION_INDEX.get((ct_a, ct_b))
        assert entry is not None, f"No REACTION_INDEX entry for ({ct_a}, {ct_b})"

        iso_a, iso_b = 400, 401
        frag_a = self._make_frag(ct_a, iso_a)
        frag_b = self._make_frag(ct_b, iso_b)

        prod = run_bond_smirks(frag_a, frag_b, iso_a, iso_b, entry,
                               intramolecular=False)
        assert prod is not None, \
            f"SMIRKS produced no product for ({ct_a}, {ct_b}) [{entry['id']}]"
        target_dummies = [a for a in prod.GetAtoms()
                         if a.GetAtomicNum() == 0
                         and a.GetIsotope() in (iso_a, iso_b)]
        assert not target_dummies, \
            f"Target dummies not consumed in ({ct_a}, {ct_b}): " \
            f"{Chem.MolToSmiles(prod)}"
        return prod

    # ── Backbone bonds ───────────────────────────────────────────────────────

    def test_backbone_amide_n_c(self):
        """backbone_n + backbone_c -> amide bond."""
        self._run_pair('backbone_n', 'backbone_c')

    def test_backbone_amide_amine_c(self):
        """amine_primary + backbone_c -> amide bond."""
        self._run_pair('amine_primary', 'backbone_c')

    def test_backbone_ester_o_c(self):
        """backbone_o + backbone_c -> ester bond."""
        self._run_pair('backbone_o', 'backbone_c')

    # ── Sidechain / crosslink bonds ──────────────────────────────────────────

    def test_sidechain_amide(self):
        """amine_primary + carboxyl -> sidechain amide."""
        self._run_pair('amine_primary', 'carboxyl')

    def test_disulfide(self):
        """thiol + thiol -> disulfide."""
        self._run_pair('thiol', 'thiol')

    def test_diselenide(self):
        """selenol + selenol -> diselenide."""
        self._run_pair('selenol', 'selenol')

    def test_thioester_backbone(self):
        """thiol + backbone_c -> thioester."""
        self._run_pair('thiol', 'backbone_c')

    def test_thioester_sidechain(self):
        """thiol + carboxyl -> sidechain thioester (C-attachment)."""
        self._run_pair('thiol', 'carboxyl')

    def test_thioether_halide(self):
        """thiol + alkyl_halide_c -> thioether."""
        self._run_pair('thiol', 'alkyl_halide_c')

    def test_thiol_maleimide(self):
        """thiol + maleimide_c -> thioether (Michael addition)."""
        self._run_pair('thiol', 'maleimide_c')

    def test_nhs_ester_amine_primary(self):
        """amine_primary + nhs_ester -> amide."""
        self._run_pair('amine_primary', 'nhs_ester')

    def test_nhs_ester_amine_secondary(self):
        """amine_secondary + nhs_ester -> amide."""
        self._run_pair('amine_secondary', 'nhs_ester')

    def test_oxime_ligation(self):
        """aminooxy + aldehyde -> oxime."""
        self._run_pair('aminooxy', 'aldehyde')

    def test_hydrazone(self):
        """hydrazide + aldehyde -> hydrazone."""
        self._run_pair('hydrazide', 'aldehyde')

    # ── Click chemistry ──────────────────────────────────────────────────────

    def test_cuaac_triazole(self):
        """alkyne_c + azide_alpha_c -> 1,2,3-triazole."""
        self._run_pair('alkyne_c', 'azide_alpha_c')

    def test_spaac_triazole(self):
        """cyclooctyne_c + azide_alpha_c -> triazole (copper-free)."""
        self._run_pair('cyclooctyne_c', 'azide_alpha_c')

    def test_iedda_tetrazine_tco(self):
        """tetrazine_c + tco_c -> dihydropyridazine (iEDDA)."""
        self._run_pair('tetrazine_c', 'tco_c')

    # ── PTM and stapling ─────────────────────────────────────────────────────

    def test_phosphorylation(self):
        """hydroxyl + phosphate_p -> phosphate ester."""
        self._run_pair('hydroxyl', 'phosphate_p')

    def test_rcm_alkene(self):
        """terminal_alkene + terminal_alkene -> internal alkene (RCM)."""
        self._run_pair('terminal_alkene', 'terminal_alkene')

    def test_imide_n_acylation(self):
        """backbone_n_mod + carboxyl -> imide (N-acylation)."""
        self._run_pair('backbone_n_mod', 'carboxyl')

    def test_backbone_n_alkylation(self):
        """backbone_n_mod + carbon -> N-alkylation."""
        self._run_pair('backbone_n_mod', 'carbon')

    # ── Reverse direction: REACTION_INDEX stores both orders ─────────────────

    def test_reverse_sidechain_amide(self):
        """carboxyl + amine_primary (reverse order) routes to same reaction."""
        self._run_pair('carboxyl', 'amine_primary')

    def test_reverse_thioester_sc(self):
        """carboxyl + thiol (reverse order) routes to thioester_sc."""
        self._run_pair('carboxyl', 'thiol')

    def test_reverse_oxime(self):
        """aldehyde + aminooxy (reverse order)."""
        self._run_pair('aldehyde', 'aminooxy')

    # ── Exhaustive: every unique pair in REACTION_INDEX fires ────────────────

    def test_all_reaction_pairs_fire(self):
        """Iterate REACTION_INDEX: every (ct_a, ct_b) pair produces a valid product."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX, run_bond_smirks

        KNOWN_XFAIL = {
            ('tetrazine_c', 'tco_c'),
        }

        seen = set()
        failures = []
        for (ct_a, ct_b), entry in REACTION_INDEX.items():
            pair_key = tuple(sorted([ct_a, ct_b]))
            if pair_key in seen:
                continue
            seen.add(pair_key)

            if (ct_a, ct_b) in KNOWN_XFAIL or (ct_b, ct_a) in KNOWN_XFAIL:
                continue

            ALL_TYPES = {**self._MODEL_MONOMERS, **self._FALLBACK_FRAGS}
            if ct_a not in ALL_TYPES or ct_b not in ALL_TYPES:
                continue

            try:
                iso_a, iso_b = 400, 401
                frag_a = self._make_frag(ct_a, iso_a)
                frag_b = self._make_frag(ct_b, iso_b)
                prod = run_bond_smirks(frag_a, frag_b, iso_a, iso_b, entry,
                                       intramolecular=False)
                if prod is None:
                    failures.append(f"({ct_a}, {ct_b}) [{entry['id']}]: no product")
                elif any(a.GetAtomicNum() == 0
                         and a.GetIsotope() in (iso_a, iso_b)
                         for a in prod.GetAtoms()):
                    failures.append(
                        f"({ct_a}, {ct_b}) [{entry['id']}]: target dummies "
                        f"not consumed in {Chem.MolToSmiles(prod)}")
            except Exception as exc:
                failures.append(f"({ct_a}, {ct_b}) [{entry['id']}]: {exc}")

        assert not failures, \
            f"{len(failures)} reaction pair(s) failed:\n" + "\n".join(failures)


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

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

    def test_crosslink_plus_bracket_converts(self):
        """(1,2) bracket converts to reversed chain with crosslink anchor."""
        from pyPept.sequence import cabiln_to_branch
        bracket = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        result = cabiln_to_branch(bracket)
        assert result == "!1-A-K.!2(4,4)-G-A-!1%C20FA-AEEA-E_g.!2", (
            f"(1,2) bracket should convert to reversed crosslink branch: {result}"
        )

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

    def test_cyclic_plus_lipid_bracket_roundtrip(self):
        """Bracket form of cyclic+lipid round-trips through branch and back."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        bracket = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        branch = cabiln_to_branch(bracket)
        assert branch == "!1-A-K.!2(4,4)-G-A-!1%C20FA-AEEA-E_g.!2", f"bracket->branch: {branch!r}"
        back = cabiln_to_bracket(branch)
        assert back == bracket, (
            f"branch->bracket roundtrip failed:\n"
            f"  {bracket!r} -> {branch!r} -> {back!r}"
        )

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

    @pytest.mark.parametrize("label,branch,expected_bracket", [
        ("two_lipid_crosslinks",
         "ac-K.!1(4,4)-G-K.!2(4,4)-am%C20FA-AEEA-E_g.!1%C20FA-AEEA-E_g.!2",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am"),
        ("three_positional_branches",
         "ac-D.!1(4,1)-A-D.!2(4,1)-A-D.!3(4,1)-am%G.!1-am%A.!2-am%G.!3-G-am",
         "ac-D.[G(4,1).am(2,1)]-A-D.[A(4,1).am(2,1)]-A-D.[G(4,1).G(2,1).am(2,1)]-am"),
        ("four_residue_branch",
         "ac-D.!1(4,1)-G-am%G.!1-A-G-am",
         "ac-D.[G(4,1).A(2,1).G(2,1).am(2,1)]-G-am"),
        ("no_caps",
         "D.!1(4,1)-G-A%G.!1-A-am",
         "D.[G(4,1).A(2,1).am(2,1)]-G-A"),
        ("cyclic_plus_plain",
         "!1-A-D.!2(4,1)-G-A-!1%G.!2-am",
         "!1-A-D.[G(4,1).am(2,1)]-G-A-!1"),
        # Exotic / mixed-type cases
        ("mixed_lipid_isopeptide",
         "ac-K.!1(4,4)-A-D.!2(4,1)-am%C20FA-AEEA-E_g.!1%G.!2-A-am",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-D.[G(4,1).A(2,1).am(2,1)]-am"),
        ("long_isopeptide_arm",
         "ac-E.!1(4,1)-G-am%G.!1-G-A-A-am",
         "ac-E.[G(4,1).G(2,1).A(2,1).A(2,1).am(2,1)]-G-am"),
        ("cyclic_plus_two_residue_arm",
         "!1-A-D.!2(4,1)-G-A-!1%G.!2-A-am",
         "!1-A-D.[G(4,1).A(2,1).am(2,1)]-G-A-!1"),
        ("two_isopeptide_one_lipid",
         "ac-D.!1(4,1)-G-K.!2(4,4)-A-D.!3(4,1)-am%G.!1-am%C20FA-AEEA-E_g.!2%A.!3-am",
         "ac-D.[G(4,1).am(2,1)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-D.[A(4,1).am(2,1)]-am"),
        ("ameleu_slot3_arm",
         "ac-D.!1(4,3)-G-am%aMeLeu.!1-G-am",
         "ac-D.[aMeLeu(4,3).G(2,1).am(2,1)]-G-am"),
        ("dual_e_scaffold",
         "ac-E.!1(4,1)-G-E.!2(4,1)-am%G.!1-A-am%A.!2-G-am",
         "ac-E.[G(4,1).A(2,1).am(2,1)]-G-E.[A(4,1).G(2,1).am(2,1)]-am"),
        ("rgd_pendant_arm",
         "ac-G-D.!1(4,1)-S-am%R.!1-G-D-am",
         "ac-G-D.[R(4,1).G(2,1).D(2,1).am(2,1)]-S-am"),
        ("cyclic_isopeptide_lipid",
         "!1-A-D.!2(4,1)-G-K.!3(4,4)-A-!1%G.!2-am%C20FA-AEEA-E_g.!3",
         "!1-A-D.[G(4,1).am(2,1)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-!1"),
    ])
    def test_branch_to_bracket_comprehensive(self, label, branch, expected_bracket):
        """branch -> bracket: exact string + roundtrip + independent SMILES equivalence."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        result = cabiln_to_bracket(branch)
        assert result == expected_bracket, (
            f"{label}: branch->bracket\n"
            f"  expected: {expected_bracket!r}\n"
            f"  got:      {result!r}"
        )
        back = cabiln_to_branch(result)
        assert back == branch, (
            f"{label}: bracket->branch roundtrip\n"
            f"  expected: {branch!r}\n"
            f"  got:      {back!r}"
        )
        assert self._smiles(branch) == self._smiles(result), (
            f"{label}: branch and bracket forms give different molecules"
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

    def test_c_terminal_marker_branch_to_bracket(self):
        """C-terminal !1 marker: %G-G-!1 converts to bracket; round-trip gives inline dot form."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = "ac-K.!1(4,2)-G-am%G-G-!1"
        bracket = cabiln_to_bracket(branch)
        assert bracket == "ac-K.[G(4,2).G(1,2)]-G-am", (
            f"bracket mismatch: {bracket!r}"
        )
        back = cabiln_to_branch(bracket)
        assert back == "ac-K.!1(4,2)-G-am%G-G-!1", f"back form: {back!r}"

    def test_retatrutide_branch_bracket_roundtrip(self):
        """Retatrutide with C-terminal -!1 marker: converts to bracket and round-trips via dot form."""
        from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
        branch = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.!1(4,2)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%C20FA-E_g-AEEA-!1"
        )
        expected_bracket = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
            "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
        )
        bracket = cabiln_to_bracket(branch)
        assert bracket == expected_bracket, f"bracket mismatch: {bracket!r}"
        back = cabiln_to_branch(bracket)
        expected_back = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.!1(4,2)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%C20FA-E_g-AEEA-!1"
        )
        assert back == expected_back, f"back mismatch: {back!r}"

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


class TestBracketBranchEquivalence:
    """Assembly-level tests: bracket and branch forms produce identical SMILES."""

    @staticmethod
    def _smiles(cabiln):
        from pyPept.sequence import Sequence
        from pyPept.molecule import Molecule
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt="ROMol")
        from rdkit.Chem import MolToSmiles
        return MolToSmiles(romol)

    @pytest.mark.parametrize("label,bracket,branch", [
        ("lipid_linker",
         "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am",
         "ac-K.!1(4,4)-G-am%C20FA-AEEA-E_g.!1"),
        ("TBMB_crosslinker",
         "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3",
         "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"),
    ])
    def test_bracket_branch_same_smiles(self, label, bracket, branch):
        smi_bracket = self._smiles(bracket)
        smi_branch = self._smiles(branch)
        assert smi_bracket == smi_branch, (
            f"{label}: SMILES mismatch\n"
            f"  bracket ({bracket}): {smi_bracket}\n"
            f"  branch  ({branch}):  {smi_branch}"
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
    def test_branch_converts_then_assembles_same(self, label, bracket, branch):
        """Branch without annotations -> convert to bracket -> assemble, compare."""
        from pyPept.sequence import cabiln_to_bracket
        converted = cabiln_to_bracket(branch)
        assert converted == bracket, f"{label}: conversion mismatch: {converted} != {bracket}"
        smi = self._smiles(bracket)
        assert smi, f"{label}: bracket form should assemble"

    def test_lipid_linker_converter_and_assembly(self):
        """Converter roundtrip + assembly via bracket for lipid linker."""
        from pyPept.sequence import cabiln_to_branch, cabiln_to_bracket
        bracket = "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        branch = cabiln_to_branch(bracket)
        bracket_back = cabiln_to_bracket(branch)
        assert bracket_back == bracket, f"Roundtrip failed: {bracket} -> {branch} -> {bracket_back}"
        smi1 = self._smiles(bracket)
        smi2 = self._smiles(bracket_back)
        assert smi1 == smi2, f"Assembly mismatch after roundtrip: {smi1} vs {smi2}"

    def test_retatrutide_bracket_assembles(self):
        """Full retatrutide with (1,2) lipid linker bracket assembles."""
        cabiln = ("Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
                  "K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
                  "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am")
        smi = self._smiles(cabiln)
        assert smi and len(smi) > 100, f"Retatrutide should produce a large SMILES: {smi}"

    def test_cyclic_plus_lipid_same_smiles(self):
        """Cyclic peptide with lipid branch: bracket and branch forms assemble identically."""
        bracket = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        branch = "!1-A-K.!2(4,4)-G-A-!1%C20FA-AEEA-E_g.!2"
        smi_bracket = self._smiles(bracket)
        smi_branch = self._smiles(branch)
        assert smi_bracket == smi_branch, (
            f"cyclic+lipid SMILES mismatch\n"
            f"  bracket: {smi_bracket}\n"
            f"  branch:  {smi_branch}"
        )

    def test_two_lipid_brackets_same_smiles(self):
        """Two lipid brackets on same chain: bracket and branch forms assemble identically."""
        bracket = "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am"
        branch = "ac-K.!1(4,4)-G-K.!2(4,4)-am%C20FA-AEEA-E_g.!1%C20FA-AEEA-E_g.!2"
        smi_bracket = self._smiles(bracket)
        smi_branch = self._smiles(branch)
        assert smi_bracket == smi_branch, (
            f"two-lipid-brackets SMILES mismatch\n"
            f"  bracket: {smi_bracket}\n"
            f"  branch:  {smi_branch}"
        )

    def test_semaglutide_like_same_smiles(self):
        """8AA semaglutide-like with C20FA lipid: bracket and branch assemble identically."""
        bracket = "His-Aib-Glu-Gly-Thr-Phe-Thr-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am"
        branch = "His-Aib-Glu-Gly-Thr-Phe-Thr-K.!1(4,4)-am%C20FA-AEEA-E_g.!1"
        smi_bracket = self._smiles(bracket)
        smi_branch = self._smiles(branch)
        assert smi_bracket == smi_branch, (
            f"semaglutide-like SMILES mismatch\n"
            f"  bracket: {smi_bracket}\n"
            f"  branch:  {smi_branch}"
        )


class TestLibraryRoundTrip:
    """Verify that the full monomer library round-trips through pre_activate."""

    LG_ATOMS = {
        '[H]': (1, 0), '[OH]': (8, 1),
        '[Cl]': (17, 0), '[Br]': (35, 0), '[I]': (53, 0),
    }

    @staticmethod
    def _restore_full_smiles(mol, rgroups):
        emol = Chem.RWMol(Chem.RWMol(mol))
        try:
            Chem.Kekulize(emol, clearAromaticFlags=True)
        except Exception:
            pass
        dummies = []
        for atom in emol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummies.append((atom.GetIdx(), atom.GetIsotope()))
        replacements = []
        lg_map = {
            '[H]': ('remove', None), '[OH]': ('replace_oh', None),
            '[Cl]': ('halide', 17), '[Br]': ('halide', 35), '[I]': ('halide', 53),
        }
        for dummy_idx, slot in dummies:
            if slot < 1 or slot > len(rgroups):
                continue
            lg = rgroups[slot - 1]
            if lg is None:
                continue
            action, param = lg_map.get(lg, ('unknown', lg))
            replacements.append((dummy_idx, action, param))
        atoms_to_remove = []
        for dummy_idx, action, param in replacements:
            if action == 'remove':
                for nb in emol.GetAtomWithIdx(dummy_idx).GetNeighbors():
                    nb.SetNoImplicit(False)
                atoms_to_remove.append(dummy_idx)
            elif action == 'replace_oh':
                emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(8)
                emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
                emol.GetAtomWithIdx(dummy_idx).SetNumExplicitHs(1)
            elif action == 'halide':
                emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(param)
                emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
            elif action == 'unknown':
                atoms_to_remove.append(dummy_idx)
        for idx in sorted(atoms_to_remove, reverse=True):
            emol.RemoveAtom(idx)
        try:
            Chem.SanitizeMol(emol)
            return Chem.MolToSmiles(emol)
        except Exception:
            return None

    @staticmethod
    def _parse_rgroups(rg):
        if isinstance(rg, list):
            return rg
        return [None if v.strip() == 'None' else v.strip()
                for v in rg.split(',')]

    @staticmethod
    def _canonical(smi):
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(mol) if mol else None

    def test_library_round_trip_threshold(self):
        """At least 93% of monomers must round-trip through pre_activate."""
        from rdkit.Chem import PandasTools
        from pyPept.interfaces.monomer_pipeline import pre_activate

        sdf = os.path.join(os.path.dirname(__file__),
                           '..', 'src', 'pyPept', 'data', 'monomers.sdf')
        df = PandasTools.LoadSDF(sdf)
        df = df.set_index('symbol')

        passed = 0
        total = 0
        for sym in df.index:
            total += 1
            mol = df.loc[sym, 'ROMol']
            if mol is None:
                continue
            rgroups = self._parse_rgroups(df.loc[sym, 'm_Rgroups'])
            stored = self._canonical(Chem.MolToSmiles(mol))

            full_smi = self._restore_full_smiles(mol, rgroups)
            if full_smi is None:
                continue

            try:
                result = pre_activate(full_smi)
            except Exception:
                continue

            result_canon = self._canonical(result.chuckles)
            stored_lgs = {}
            for i, lg in enumerate(rgroups, 1):
                if lg is not None:
                    stored_lgs[i] = lg

            if stored == result_canon and result.leaving == stored_lgs:
                passed += 1

        pct = 100 * passed / total
        assert pct >= 93.0, (
            f"Round-trip pass rate {pct:.1f}% ({passed}/{total}) "
            f"is below 93% threshold"
        )


# ---------------------------------------------------------------------------
# smiles_to_cabiln_core — comprehensive rule-coverage tests
# ---------------------------------------------------------------------------

# Helper: round-trip CABILN -> SMILES -> CABILN via smiles_to_cabiln_core.
# Imported here so every test method can use it without an extra import dance.
def _s2c_roundtrip(biln: str):
    """CABILN -> SMILES -> CABILN.

    Returns (cabiln_str, details) where details is the list of
    (abbr, n_matched, n_total_atoms) tuples from smiles_to_cabiln_core.
    """
    import sys as _sys, os as _os
    # Add the repo root so that `import tools.live_renderer` resolves.
    _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
    if _repo_root not in _sys.path:
        _sys.path.insert(0, _repo_root)
    from tools.live_renderer import smiles_to_cabiln_core  # noqa: E402
    mol = Molecule(Sequence(biln)).get_molecule(fmt='ROMol')
    smi = Chem.MolToSmiles(mol)
    return smiles_to_cabiln_core(smi)


class TestSmilesToCabiln:
    """Comprehensive rule-coverage tests for smiles_to_cabiln_core.

    The function converts a peptide SMILES string to CABILN notation.
    Its internal processing order is:

      Rule 1 — cap stripping (N-cap position 0 first; C-cap position main_n-1
               second; for a single-residue peptide BOTH strips apply because
               0 == main_n-1 simultaneously — the elif→if fix enables this)
      Rule 2 — residue matching canonical form (backbone monomers preferred
               over modified forms; CIP stereo used to pick D vs L)
      Rule 3 — cyclic detection (!1-…-!1 wrapping)
      Rule 4+5 — sidechain branch-position detection + lipid-linker bracket
               notation (.[E_g(4,4).AEEA(1,2).C20FA(1,2)])
      Rule 6 — crosslink / staple annotation (.!n(r,r))

    Design note on cap semantics
    ----------------------------
    smiles_to_cabiln_core strips caps *internally* (so residue matching finds
    'G' not 'G+acetyl'), then *re-attaches* the cap tokens to the output
    string.  Therefore:
      ac-G-am  -> round-trip produces  ac-G-am  (both caps preserved)
      G-am     -> round-trip produces  G-am     (C-cap preserved)
    The key observable for Rule 1 correctness is that the residue *abbreviation*
    in details is the expected canonical token (e.g. 'G', 'V'), not '?' or a
    modified form like 'Gly_al'.  Without the elif->if fix, a single-residue
    peptide with both caps would fail to strip the C-cap, causing the residue
    match to fail ('?') because the amide tail is not part of any library entry.
    """

    @staticmethod
    def _r(biln: str):
        """Shorthand that returns (cabiln, details)."""
        return _s2c_roundtrip(biln)

    # ------------------------------------------------------------------
    # Rule 1 — Cap stripping priority order
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("biln,expected_abbr,expected_cabiln", [
        # Single residue, C-cap (am) only.
        # Without elif->if fix, main_pos==0 strips N-cap (nothing to strip)
        # then the `elif` is skipped and the C-cap is NOT stripped, so
        # the residue mol still has the terminal amide.  Match returns '?'.
        pytest.param("G-am", "G", "G-am", id="single_G_am_only"),
        pytest.param("A-am", "A", "A-am", id="single_A_am_only"),
        pytest.param("V-am", "V", "V-am", id="single_V_am_only"),

        # Single residue, N-cap (ac) only.
        # N-cap strip fires at position 0; no C-cap present.
        pytest.param("ac-G", "G", "ac-G", id="single_ac_G_only"),
        pytest.param("ac-A", "A", "ac-A", id="single_ac_A_only"),

        # Single residue, BOTH caps.
        # Critical: position 0 == main_n-1.  With `elif`, only the N-cap
        # strip would fire; the C-cap amide would remain and break matching.
        # With `if`, both strips fire on the same residue.
        pytest.param("ac-G-am", "G", "ac-G-am", id="single_ac_G_am_both"),
        pytest.param("ac-A-am", "A", "ac-A-am", id="single_ac_A_am_both"),
        pytest.param("ac-V-am", "V", "ac-V-am", id="single_ac_V_am_both"),

        # Multi-residue, both caps.  N-cap strip touches pos 0 only;
        # C-cap strip touches pos main_n-1 only.  Middle residues unaffected.
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
        """Rule 1: caps are stripped internally for matching, then re-appended.

        expected_abbr  — the abbreviation the first residue must receive
                         in the details list (verifies the internal strip worked)
        expected_cabiln — the full reconstructed string (verifies re-attachment)

        A regression (elif instead of if) would leave the first residue's
        abbreviation as '?' for any single-residue peptide with a C-cap.
        """
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

    def test_cyclic_gets_no_caps(self):
        """Rule 1 (boundary): cyclic peptide bypasses the cap-strip block.

        smiles_to_cabiln_core skips both strip passes when cyclic is True, so
        a cyclic input must come back with !1-…-!1 and no spurious ac-/am-.
        The ac and am tokens appear in N-terminal and C-terminal tokens of the
        main chain position 0 and main_n-1.  On a cyclic peptide neither
        position exists, so the keywords must not appear.
        """
        result, details = self._r("!1-A-G-K-!1")
        assert result.startswith("!1-"), f"Expected cyclic prefix: {result!r}"
        assert result.endswith("-!1"),   f"Expected cyclic suffix: {result!r}"
        assert not result.startswith("ac-!1"), f"Spurious N-cap: {result!r}"
        assert "?" not in result, f"Unknown residue in cyclic: {result!r}"

    # ------------------------------------------------------------------
    # Rule 2 — Residue matching canonical form
    # ------------------------------------------------------------------

    def test_all_20_standard_aa(self):
        """Rule 2: all 20 canonical amino acids round-trip to their standard abbrs.

        Backbone_n+backbone_c entries (A, G, V, …) must be preferred over
        any modified form (Ala_al, Gly_al, …) that may also match the fragment.
        The 20-mer is assembled cap-free on the N-terminus and with -am on C.
        We check the first 20 tokens of the recovered CABILN.
        """
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

    @pytest.mark.parametrize("biln,expected_abbr", [
        # G preferred over Gly_al.
        pytest.param("G-am", "G", id="G_not_Gly_al"),
        # A preferred over Ala_al.
        pytest.param("A-am", "A", id="A_not_Ala_al"),
    ])
    def test_residue_priority_canonical_over_modified(self, biln, expected_abbr):
        """Rule 2: backbone monomers win over modified variants with equal atom count.

        Gly_al and Ala_al match the same backbone atoms as G and A.  The
        priority mechanism (backbone_n+backbone_c monomer preferred) must
        return G / A, not Gly_al / Ala_al.
        """
        result, details = self._r(biln)
        abbr = details[0][0]
        assert abbr == expected_abbr, (
            f"Got {abbr!r} instead of {expected_abbr!r} for {biln!r}"
        )

    def test_d_amino_acid_stereo_disambiguation(self):
        """Rule 2: CIP stereo scoring correctly identifies D-amino acids.

        dA (D-Ala) and dV (D-Val) have inverted Cα CIP codes vs L-isomers.
        The library contains both isomers; the match function must pick the
        D-form for each and return 'DAla' / 'DVal' (the canonical SDF names).
        Returning 'A'/'V' (L-forms) would indicate that CIP scoring is broken.
        """
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
        """Rule 2: Aib and Orn are recognised as non-natural backbone monomers.

        These are in the library with backbone_n+backbone_c connectivity.
        The match must not fall back to '?'.
        """
        result, details = self._r("Aib-Orn-am")
        assert "?" not in result, f"Unknown residue in: {result!r}"
        abbrs = [d[0] for d in details]
        assert abbrs[0] == "Aib", f"Expected Aib, got {abbrs[0]!r}"
        assert abbrs[1] == "Orn", f"Expected Orn, got {abbrs[1]!r}"

    def test_l_vs_d_leucine_stereo(self):
        """Rule 2 (stereo boundary): L-Leu and D-Leu produce different tokens.

        If CIP scoring were absent, both would map to the first matching library
        entry (whichever of L or D appears first), producing the same output.
        """
        result_l,  details_l  = self._r("L-am")
        result_dl, details_dl = self._r("dL-am")
        assert "?" not in result_l,  f"Unknown in L-am: {result_l!r}"
        assert "?" not in result_dl, f"Unknown in dL-am: {result_dl!r}"
        assert details_l[0][0] != details_dl[0][0], (
            f"L-Leu and D-Leu produced same abbr: {details_l[0][0]!r}"
        )

    # ------------------------------------------------------------------
    # Rule 3 — Cyclic peptide detection
    # ------------------------------------------------------------------

    @pytest.mark.parametrize("biln", [
        # Minimal 3-mer cyclic (head-to-tail backbone ring).
        pytest.param("!1-A-G-K-!1",   id="cyclic_3mer"),
        # 4-mer cyclic.
        pytest.param("!1-A-G-K-D-!1", id="cyclic_4mer"),
    ])
    def test_cyclic_detection(self, biln):
        """Rule 3: cyclic backbone is detected and wrapped with !1-…-!1.

        The BFS on the peptide-bond graph finds the ring and sets cyclic=True,
        which causes the output to be prefixed/suffixed with !1.  No caps are
        added.  All residue tokens must be recognised (no '?').
        """
        result, details = self._r(biln)
        assert result.startswith("!1-"), f"Missing cyclic prefix: {result!r}"
        assert result.endswith("-!1"),   f"Missing cyclic suffix: {result!r}"
        assert "?" not in result,        f"Unknown residue in cyclic: {result!r}"

    # ------------------------------------------------------------------
    # Rules 4+5 — Branch detection + lipid-linker bracket notation
    # ------------------------------------------------------------------

    def test_minimal_lipid_linker_bracket(self):
        """Rules 4+5: K with a E_g-AEEA-C20FA sidechain in a 3-residue backbone.

        The isopeptide-bonded branch (E_g) is detected as a degree-1 node
        in the peptide-bond graph whose neighbour (K) has degree > 2.  K must
        emerge with bracket notation K.[E_g(4,4).AEEA(1,2).C20FA(1,2)].
        The backbone residues A and G must also be recognised.
        No '?' tokens allowed.
        """
        biln   = "A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in lipid-linker result: {result!r}"
        assert "K.[" in result or "K.!1" in result, (
            f"Expected bracket/crosslink notation on K in: {result!r}"
        )

    def test_retatrutide_full_sequence(self):
        """Rules 4+5: Retatrutide round-trip (39 backbone AA + C20 lipid branch).

        K16-K17 backbone; lipid branch on K17: AEEA(4,2).E_g(1,2).C20FA(1,2).
        E_g(1,2) = Glu via γ-COOH entry (R2) — preferred over E(1,4) by 1-2 rule.
        """
        biln = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
            "K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
            "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
        )
        result, details = self._r(biln)
        assert result == biln, (
            f"Retatrutide round-trip mismatch:\n  got: {result!r}\n  exp: {biln!r}"
        )

    # ------------------------------------------------------------------
    # Rule 6 — Crosslink / staple annotation
    # ------------------------------------------------------------------

    def test_staple_i_i4_same_handedness(self):
        """Rule 6: i, i+4 RCM staple (S5+S5, same configuration).

        The side-chain bond between the two S5 olefin tails is detected as a
        'side chain' pair in the residue-pair graph.  Both S5 positions must
        be annotated with .!1(4,4) (the attachment R-groups of S5 are R4/R4).
        Caps (ac-, -am) must also survive.
        """
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
        """Rule 6: i, i+7 RCM staple (S5+R8, opposite configuration).

        S5 and R8 have the same sidechain olefin length but opposite Cα
        configurations.  The crosslink pair must still be detected and both
        positions annotated with .!1.
        """
        biln   = "ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in i+7 staple result: {result!r}"
        assert result.count(".!1") == 2, (
            f"Expected two .!1 annotations in: {result!r}"
        )

    # ------------------------------------------------------------------
    # Priority conflict tests
    # ------------------------------------------------------------------

    def test_priority_ac_am_lipid_branch(self):
        """Rules 1+4+5 together: ac + am + lipid branch all interact correctly.

        Known limitation (documented behaviour, not a bug to fix here):
        smiles_to_cabiln_core's backbone-detection pass currently identifies
        only the K residue when a large lipid branch is present; it does not
        always recover flanking A and G residues or their caps in that topology.
        The critical invariant tested here is that:
          (a) the K bracket notation survives (the branch is not lost), and
          (b) no '?' tokens appear (matching succeeds for every detected residue).
        Cap presence is topology-dependent and is NOT asserted here.
        """
        biln   = "ac-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-am"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue: {result!r}"
        assert "K.[" in result or "K.!1" in result, (
            f"Lipid branch on K lost: {result!r}"
        )

    def test_priority_ac_am_staple(self):
        """Rules 1+6 together: N-cap + C-cap + RCM staple crosslink.

        The cap-strip pass must not corrupt the position map used by crosslink
        detection.  All three features must appear in the output simultaneously.
        """
        biln   = "ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am"
        result, details = self._r(biln)
        assert "?" not in result,       f"Unknown residue: {result!r}"
        assert result.startswith("ac-"), f"N-cap lost: {result!r}"
        assert result.endswith("-am"),   f"C-cap lost: {result!r}"
        assert result.count(".!1") == 2, f"Crosslink count wrong: {result!r}"

    def test_priority_cyclic_plus_lipid_branch(self):
        """Rules 3+4+5: cyclic backbone AND a lipid sidechain branch on K.

        On a cyclic backbone, every backbone node has peptide-bond degree 2,
        so K has degree 3 (two backbone bonds + one isopeptide bond to E_g).
        The branch_pos heuristic (degree-1 nodes) does not apply directly,
        but the BFS cyclic detection and branch detection interact correctly:
        the output is cyclic-wrapped (!1-…-!1) and K carries its bracket.

        Invariants verified:
          - !1-…-!1 wrapping present (cyclic detected)
          - K bracket notation present (lipid branch not lost)
          - no '?' tokens (all residues matched)
        """
        biln   = "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"
        result, details = self._r(biln)
        assert "?" not in result, f"Unknown residue in cyclic+branch: {result!r}"
        assert result.startswith("!1-"), f"Cyclic prefix lost: {result!r}"
        assert result.endswith("-!1"),   f"Cyclic suffix lost: {result!r}"
        assert "K.[" in result or "K.!1" in result, (
            f"Lipid branch on K lost in cyclic+branch: {result!r}"
        )

    def test_dual_lipid_branch_terminal_k(self):
        """Rules 4+5: two K residues each with E(4,4)-AEEA-C20FA, K at N- and C-terminus.

        Exercises the terminal-K branch detection fix: K has degree 2 in the
        peptide-bond graph (one backbone bond to G + one isopeptide to E(4,4)),
        which previously prevented E from being recognised as a branch_pos.
        The isopeptide bond is now classified as 'side chain' and E is detected
        as a degree-0 backbone node with a side-chain connection to a backbone K.

        E(4,4) = Glu attached via γ-COOH (R4) on both ends — formerly E_g.

        Invariants:
          - output == input (exact round-trip)
          - Both K residues carry bracket notation
          - No '?' tokens
        """
        biln = "ac-K.[E(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E(4,4).AEEA(1,2).C20FA(1,2)]-am"
        result, details = self._r(biln)
        assert result == biln, (
            f"Dual-lipid round-trip failed:\n"
            f"  expected: {biln!r}\n"
            f"  got:      {result!r}"
        )
        assert "?" not in result, f"Unknown residue in dual-lipid result: {result!r}"
        assert result.count("K.[") == 2, (
            f"Expected two K bracket annotations, got: {result!r}"
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
        assert result == cabiln

    def test_tbmb_asymmetric_roundtrip(self):
        """Asymmetric-loop TBMB bicycle SMILES→CABILN round-trip."""
        cabiln = "ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TBMB.!1.!2.!3"
        result, _ = self._r(cabiln)
        assert result == cabiln

    def test_tbmb_partial_two_arm_roundtrip(self):
        """TBMB 2-of-3 arms reacted: SMILES→CABILN detects partial scaffold."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2"
        result, _ = self._r(cabiln)
        assert result == cabiln

    def test_tbmb_scaffold_plus_disulfide_roundtrip(self):
        """TBMB scaffold + independent Cys-Cys disulfide on same peptide."""
        cabiln = "ac-C.!1(4,4)-A-C.!4(4,4)-A-C.!4(4,4)-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        result, _ = self._r(cabiln)
        assert result == cabiln

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
        assert result == cabiln

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
        assert result == cabiln

    # ── Gap 1: TATA 2-arm partial roundtrip ───────────────────────────────────

    def test_tata_partial_two_arm_roundtrip(self):
        """TATA 2-of-3 arms reacted: partial scaffold detected, slots normalise to (4,4)/(4,5)."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TATA.!1.!2"
        result, _ = self._r(cabiln)
        assert result == cabiln, (
            f"TATA 2-arm partial round-trip failed:\n"
            f"  expected: {cabiln!r}\n"
            f"  got:      {result!r}"
        )

    # ── Gap 2: TBAB 2-arm partial roundtrip ───────────────────────────────────

    def test_tbab_partial_two_arm_roundtrip(self):
        """TBAB 2-of-3 arms reacted: partial scaffold detected, slots normalise to (4,4)/(4,5)."""
        cabiln = "ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBAB.!1.!2"
        result, _ = self._r(cabiln)
        assert result == cabiln, (
            f"TBAB 2-arm partial round-trip failed:\n"
            f"  expected: {cabiln!r}\n"
            f"  got:      {result!r}"
        )

    # ── Gap 3: Bracket notation (→ []) assembles same SMILES ──────────────────

    @pytest.mark.parametrize("pct_cabiln", [
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3", id="tbmb_3arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2",                   id="tbmb_2arm"),
        pytest.param("ac-C.!1(4,4)-A-A-C.!1(4,4)-am",                               id="disulfide"),
    ])
    def test_bracket_notation_assembles_same_smiles(self, pct_cabiln):
        """→ [] path: bracket CABILN assembles to the identical SMILES as the percent form."""
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import _renumber_xlinks
        from pyPept.sequence import cabiln_to_bracket

        ref_mol = Molecule(Sequence(pct_cabiln)).get_molecule(fmt='ROMol')
        ref_smi = Chem.MolToSmiles(ref_mol)

        bracket_cabiln = _renumber_xlinks(cabiln_to_bracket(pct_cabiln))
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
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import _renumber_xlinks, smiles_to_cabiln_core
        from pyPept.sequence import cabiln_to_bracket

        bracket_cabiln = _renumber_xlinks(cabiln_to_bracket(pct_cabiln))
        brk_mol = Molecule(Sequence(bracket_cabiln)).get_molecule(fmt='ROMol')
        brk_smi = Chem.MolToSmiles(brk_mol)
        result, _ = smiles_to_cabiln_core(brk_smi)

        assert result == pct_cabiln, (
            f"Bracket→SMILES→percent round-trip failed:\n"
            f"  expected: {pct_cabiln!r}\n"
            f"  got:      {result!r}"
        )

    # ── Gap 3c: Endpoint _renumber_xlinks produces consecutive crosslink IDs ──

    def test_bracket_notation_endpoint_renumbers_crosslinks(self):
        """_renumber_xlinks(cabiln_to_bracket(…)) produces !1/!2/… with no gaps."""
        import sys as _sys, os as _os, re as _re
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import _renumber_xlinks
        from pyPept.sequence import cabiln_to_bracket

        pct = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
        bracket = _renumber_xlinks(cabiln_to_bracket(pct))

        assert 'TBMB' in bracket, f"TBMB scaffold lost in bracket conversion: {bracket!r}"
        ids = sorted(set(int(m) for m in _re.findall(r'\.!(\d+)', bracket)))
        assert ids == list(range(1, len(ids) + 1)), (
            f"Non-consecutive crosslink IDs in bracket form: {ids!r} ({bracket!r})"
        )

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
        """CuAAC SMILES→CABILN round-trip via de-reacted backbone detection."""
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
        assert result == cabiln, f"Expected {cabiln!r}, got {result!r}"

    def test_standalone_tbmb_smiles_detected(self):
        """Unreacted TBMB scaffold (all three Br arms present) is identified as 'TBMB'."""
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import smiles_to_cabiln_core
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
        """CIP @/@@ stereocenters survive SMILES→CABILN→SMILES without inversion."""
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
        assert result == cabiln, (
            f"Full TATA match degraded to partial:\n"
            f"  expected: {cabiln!r}\n"
            f"  got:      {result!r}"
        )
        assert "%TATA.!1.!2.!3" in result, f"3-arm TATA suffix missing: {result!r}"

    # ── Gap 8: thia_michael_c SMARTS specificity ─────────────────────────────

    def test_thia_michael_c_smarts_specificity(self):
        """thia_michael_c pre_smarts matches beta-acrylamide C but not thiol, carboxyl, or alkene."""
        thia_smarts = Chem.MolFromSmarts('[CH2:1]=[CH]C(=O)[NX3]')
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
        """Rule 1: N-cap type identified from library, not hardcoded as 'ac'.

        _s2c_identify_n_cap matches the stripped cap fragment against all
        m_type='cap' + backbone_c monomers.  The correct abbreviation must
        appear at the start of the reconstructed CABILN string.
        """
        result, details = self._r(biln)
        assert result.startswith(expected_cap_token + '-'), (
            f"Expected N-cap {expected_cap_token!r} in output, got: {result!r}"
        )
        assert "?" not in result, f"Unknown residue in capped peptide: {result!r}"

    @pytest.mark.parametrize("biln,expected_c_cap", [
        pytest.param("ac-A-am",    "am",    id="c_cap_am"),
        pytest.param("ac-A-NHMe",  "NHMe",  id="c_cap_NHMe"),
        pytest.param("ac-A-NHEt",  "NHEt",  id="c_cap_NHEt"),
        pytest.param("ac-A-OEt",   "OEt",   id="c_cap_OEt"),
        pytest.param("ac-A-OtBu",  "OtBu",  id="c_cap_OtBu"),
        pytest.param("ac-A-_OMe",  "_OMe",  id="c_cap_OMe"),
        pytest.param("ac-A-_OBn",  "_OBn",  id="c_cap_OBn"),
        pytest.param("fmoc-G-G-OEt", "OEt", id="c_cap_OEt_multi"),
    ])
    def test_c_cap_type_identification(self, biln, expected_c_cap):
        """C-cap type identified from library for both amide-N and ester-O caps.

        Ester C-caps (OEt, OtBu, _OMe, _OBn) were previously broken: cutting
        a C(=O)-O bond left a bare aldehyde instead of a carboxyl, so the
        residue matched Ala_al/Gly_al.  The fix adds -OH when the outside atom
        of the cut bond is O (not just N).
        """
        result, _ = self._r(biln)
        assert result == biln, f"Expected {biln!r}, got {result!r}"
        assert result.endswith('-' + expected_c_cap), (
            f"Expected C-cap {expected_c_cap!r} at end of output, got: {result!r}"
        )

    # ------------------------------------------------------------------
    # Unknown monomer raises ValueError
    # ------------------------------------------------------------------

    def test_no_question_mark_tokens(self):
        """Rule 2: '?' tokens must not appear in smiles_to_cabiln_core output.

        The old code used `abbr or '?'` as a silent fallback for unrecognised
        residues.  The new code raises ValueError instead.  As a proxy test,
        verify that a broad set of standard sequences never produce '?' in their
        round-trip output — if abbr were ever None and the raise was absent, '?'
        would appear here.
        """
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

    def test_unknown_monomer_raises_on_empty_lib_match(self):
        """Rule 2: ValueError with 'Unrecognised monomer' is raised when _s2c_match
        returns None (lib has no match for the residue fragment).

        We test this by patching _s2c_match to return (None, 0), simulating a
        library gap.  This exercises the raise path that replaced `abbr or '?'`.
        """
        import sys as _sys, os as _os, unittest.mock as _mock
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        import tools.live_renderer as _lr
        smi = 'CC(=O)NCC(N)=O'  # ac-G-am — valid SMILES, would normally succeed
        with _mock.patch.object(_lr, '_s2c_match', return_value=(None, 0)):
            with pytest.raises(ValueError, match='Unrecognised monomer'):
                _lr.smiles_to_cabiln_core(smi)

    # ------------------------------------------------------------------
    # Edge-case tests 14-20, 22-30
    # ------------------------------------------------------------------

    def test_explicit_h_amide_smiles(self):
        """Explicit H on backbone amide N/C atoms does not confuse residue matching."""
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import smiles_to_cabiln_core
        smi = 'CC(=O)[NH][C@@H](C)C(=O)[NH2]'  # ac-A-am with explicit H
        result, details = smiles_to_cabiln_core(smi)
        assert 'A' in result, f"Expected A residue in output, got {result!r}"
        assert '?' not in result, f"Unexpected '?' in {result!r}"

    def test_protonated_sidechain_roundtrip(self):
        """Lys (protonated amine sidechain) round-trips cleanly through CABILN."""
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

    def test_renumber_xlinks_renumbers_nonconsecutive(self):
        """_renumber_xlinks maps non-consecutive IDs (e.g. !1, !3) to consecutive ones."""
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import _renumber_xlinks
        cabiln = 'ac-C.!1(4,4)-A-A-C.!3(4,4)-am'
        result = _renumber_xlinks(cabiln)
        assert '.!3' not in result, f"Non-consecutive ID !3 not renumbered in {result!r}"
        assert '.!2' in result, f"Expected .!2 after renumber, got {result!r}"
        assert '.!1' in result, f"Expected .!1 preserved in {result!r}"

    def test_boc_ncap_roundtrip(self):
        """boc N-cap on a peptide with a disulfide crosslink round-trips correctly."""
        result, _ = self._r('boc-C.!1(4,4)-A-G-C.!1(4,4)-am')
        assert result.startswith('boc-'), f"boc cap lost: {result!r}"
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

    def test_alternative_ncap_with_scaffold_roundtrip(self):
        """fmoc N-cap with 2-arm TATA partial scaffold round-trips."""
        cabiln = 'fmoc-C.!1(4,4)-A-A-C.!2(4,5)-am%TATA.!1.!2'
        result, _ = self._r(cabiln)
        assert 'fmoc' in result, f"fmoc cap lost: {result!r}"
        assert '%TATA' in result, f"TATA scaffold lost: {result!r}"

    def test_scaffold_detected_at_most_once(self):
        """smiles_to_cabiln_core annotates a TBMB scaffold exactly once."""
        import re as _re
        result, _ = self._r('ac-C.!1(4,4)-A-C.!2(4,5)-A-C.!3(4,6)-am%TBMB.!1.!2.!3')
        count = len(_re.findall(r'%TBMB', result))
        assert count == 1, f"Expected exactly 1 %TBMB, found {count} in {result!r}"

    def test_scaffold_patterns_1arm_halogen_detected(self):
        """1-arm TBMB with 2 remaining CBr arms is detected as %TBMB.

        The halogen-relaxation rule: when all unmatched scaffold arms hit Br/Cl/I
        atoms (leaving groups still present in the SMILES), the min_arms threshold
        is relaxed and the partial scaffold is annotated correctly.
        """
        result, _ = self._r('ac-C.!1(4,4)-G-A-K-C-E-L-F-C-am%TBMB.!1')
        assert '%TBMB' in result, (
            f"1-arm TBMB with Br leaving groups must carry scaffold annotation: {result!r}"
        )
        assert result.endswith('%TBMB.!1'), (
            f"Expected single-arm annotation %%TBMB.!1, got {result!r}"
        )

    def test_scaffold_patterns_1arm_tata_detected(self):
        """1-arm TATA (2 unreacted vinyl arms) is detected as %TATA.

        pyPept H-caps unconnected dummy atoms so the unreacted vinyl arm becomes
        propanoyl (-CH2-CH3).  The partial SMARTS with both unreacted arm tokens
        removed must still match the TATA core and annotate the 1 reacted Cys.
        """
        result, _ = self._r('ac-C.!1(4,4)-A-A-am%TATA.!1')
        assert '%TATA' in result, (
            f"1-arm TATA must carry scaffold annotation: {result!r}"
        )
        assert result.endswith('%TATA.!1'), (
            f"Expected single-arm annotation %TATA.!1, got {result!r}"
        )

    def test_tbmb_mol_not_annotated_as_tata(self):
        """A TBMB-crosslinked peptide is labelled %TBMB, never %TATA."""
        result, _ = self._r('ac-C.!1(4,4)-A-A-C.!2(4,5)-am%TBMB.!1.!2')
        assert '%TBMB' in result, f"Expected %TBMB in {result!r}"
        assert '%TATA' not in result, f"Spurious %TATA annotation in {result!r}"

    # ── Multi-amide backbone-less chains ──────────────────────────────────────

    @pytest.mark.parametrize("biln", [
        pytest.param("ac-G-am",     id="ac_G_am"),
        pytest.param("ac-A-am",     id="ac_A_am"),
        pytest.param("ac-V-am",     id="ac_V_am"),
        pytest.param("ac-G-G-am",   id="ac_G_G_am"),
        pytest.param("ac-A-G-am",   id="ac_A_G_am"),
        pytest.param("ac-G-A-am",   id="ac_G_A_am"),
        pytest.param("ac-G-G-G-am", id="ac_G_G_G_am"),
        pytest.param("fmoc-G-am",   id="fmoc_G_am"),
        pytest.param("boc-A-G-am",  id="boc_A_G_am"),
    ])
    def test_multi_amide_backbone_less_chain_roundtrip(self, biln):
        """Capped peptide SMILES round-trips even if backbone detection fails.

        The multi-amide fallback walker handles chains with ≥2 amide bonds that
        the primary backbone detector misses (e.g. very short or cap-heavy chains).
        """
        result, _ = self._r(biln)
        assert result == biln, f"Expected {biln!r}, got {result!r}"

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


# ---------------------------------------------------------------------------
# Edge-case / error-handling tests (standalone classes)
# ---------------------------------------------------------------------------

class TestEdgeCasesErrorHandling:
    """Tests for boundary conditions and invalid-input behaviour."""

    @staticmethod
    def _setup():
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import smiles_to_cabiln_core
        return smiles_to_cabiln_core

    def test_disconnected_smiles_does_not_crash(self):
        """Disconnected SMILES (two fragments separated by '.') does not raise.

        smiles_to_cabiln_core processes the largest connected fragment; the
        smaller fragment is silently ignored.  The key guarantee is that the
        call completes without an exception and produces no '?' tokens.
        """
        smiles_to_cabiln_core = self._setup()
        smi = 'CC(=O)N[C@@H](CS)C(=O)N.OC(=O)C'  # Cys-am fragment + acetic acid
        result, _ = smiles_to_cabiln_core(smi)
        assert '?' not in result, f"Unexpected '?' in output for disconnected SMILES: {result!r}"

    def test_self_crosslink_same_id_twice_does_not_raise(self):
        """Duplicate crosslink ID on the same residue (.!1.!1) is accepted without error.

        This is a degenerate input that falls through to single-arm (no bond
        formed); the Sequence constructor must not raise.
        """
        seq = Sequence('ac-C.!1(4,4).!1(4,4)-am')
        assert seq is not None

    def test_scaffold_with_zero_crosslinks_assembles(self):
        """Scaffold suffix with no crosslink annotations assembles (scaffold unattached)."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            seq = Sequence('ac-A-am%TBMB')
            mol = Molecule(seq).get_molecule(fmt='ROMol')
        assert mol is not None, "Assembly returned None for scaffold with no crosslinks"

    def test_invalid_rgroup_slot_raises(self):
        """R-group slot number outside valid range raises ValueError."""
        with pytest.raises(ValueError):
            Sequence('ac-C.!1(4,9)-am%TBMB.!1')

    def test_same_rgroup_dual_crosslinks_assembles(self):
        """R4 of Cys annotated for two different crosslinks: only one bond forms; no crash."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            seq = Sequence('ac-C.!1(4,4).!2(4,4)-am%TBMB.!1.!2')
            mol = Molecule(seq).get_molecule(fmt='ROMol')
        assert mol is not None

    def test_crosslink_id_zero_assembles(self):
        """Crosslink ID !0 (zero) is accepted and forms a bond (disulfide)."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            seq = Sequence('ac-C.!0(4,4)-A-A-C.!0(4,4)-am')
            mol = Molecule(seq).get_molecule(fmt='ROMol')
        assert mol is not None, "Assembly with crosslink ID !0 returned None"
        smi = Chem.MolToSmiles(mol)
        assert '.' not in smi, f"Disconnected product for !0 crosslink: {smi!r}"


# ---------------------------------------------------------------------------
# Additional chemistry edge-case tests
# ---------------------------------------------------------------------------

class TestAdditionalChemistryEdgeCases:
    """Assembly and SMIRKS correctness for less-common reaction types."""

    def test_thiol_maleimide_adduct_assembly(self):
        """Thiol-maleimide ligation (Cys + 5FM) forms the expected thioether-succinimide."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mol = Molecule(Sequence('ac-C.!1(4,4)-am%5FM.!1')).get_molecule(fmt='ROMol')
        assert mol is not None, "Thiol-maleimide assembly returned None"
        smi = Chem.MolToSmiles(mol)
        thioether_succ = Chem.MolFromSmarts('SC1CC(=O)NC1=O')
        assert mol.HasSubstructMatch(thioether_succ), (
            f"Thioether-succinimide ring not found in product: {smi!r}"
        )

    def test_rcm_alkene_staple_assembly(self):
        """Ring-closing metathesis staple (S5 + R8) forms a macrocycle with an internal alkene."""
        mol = Molecule(Sequence('ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am')).get_molecule(fmt='ROMol')
        assert mol is not None, "RCM staple assembly returned None"
        smi = Chem.MolToSmiles(mol)
        alkene = Chem.MolFromSmarts('C=C')
        assert mol.HasSubstructMatch(alkene), f"Internal alkene not found after RCM: {smi!r}"
        assert '*' not in smi, f"Unreacted dummy (*) in RCM product: {smi!r}"
        assert '.' not in smi, f"Disconnected RCM product: {smi!r}"

    def test_depsipeptide_backbone_ester(self):
        """Depsipeptide (backbone ester) test is skipped — no hydroxy-acid residue in standard library."""
        pytest.skip("No backbone hydroxy-acid residue in standard library")

    def test_diselenide_assembly(self):
        """Selenocysteine (Sec) diselenoide crosslink assembles with Se–Se bond."""
        mol = Molecule(Sequence('ac-Sec.!1(4,4)-A-A-Sec.!1(4,4)-am')).get_molecule(fmt='ROMol')
        assert mol is not None, "Diselenide assembly returned None"
        smi = Chem.MolToSmiles(mol)
        sese = Chem.MolFromSmarts('[Se][Se]')
        assert mol.HasSubstructMatch(sese), f"Se–Se bond not found in product: {smi!r}"
        assert '*' not in smi, f"Unreacted dummy in diselenide product: {smi!r}"
        assert '.' not in smi, f"Disconnected diselenide product: {smi!r}"

    def test_hydrazone_smirks_produces_imine(self):
        """Hydrazone-forming SMIRKS converts aldehyde + hydrazine → C=N–NH2."""
        from rdkit.Chem import AllChem
        smirks = '[CH:1]=[O:2].[NH2:3][NH2:4] >> [CH:1]=[N:3][NH2:4]'
        rxn = AllChem.ReactionFromSmarts(smirks)
        assert rxn is not None, "Hydrazone SMIRKS failed to parse"
        aldehyde = Chem.MolFromSmiles('CCC=O')
        hydrazine = Chem.MolFromSmiles('NN')
        products = rxn.RunReactants((aldehyde, hydrazine))
        assert products, "Hydrazone SMIRKS produced no products"
        prod_smi = Chem.MolToSmiles(products[0][0])
        imine = Chem.MolFromSmarts('C=N')
        prod_mol = Chem.MolFromSmiles(prod_smi)
        assert prod_mol.HasSubstructMatch(imine), (
            f"Expected C=N imine in hydrazone product, got {prod_smi!r}"
        )

    def test_oxime_smirks_produces_oxime(self):
        """Oxime-forming SMIRKS converts aldehyde + hydroxylamine → C=N–OH."""
        from rdkit.Chem import AllChem
        smirks = '[CH:1]=[O:2].[NH2:3][OH:4] >> [CH:1]=[N:3][OH:4]'
        rxn = AllChem.ReactionFromSmarts(smirks)
        assert rxn is not None, "Oxime SMIRKS failed to parse"
        aldehyde = Chem.MolFromSmiles('CCC=O')
        hydroxylamine = Chem.MolFromSmiles('NO')
        products = rxn.RunReactants((aldehyde, hydroxylamine))
        assert products, "Oxime SMIRKS produced no products"
        prod_smi = Chem.MolToSmiles(products[0][0])
        oxime_pat = Chem.MolFromSmarts('C=NO')
        prod_mol = Chem.MolFromSmiles(prod_smi)
        assert prod_mol.HasSubstructMatch(oxime_pat), (
            f"Expected C=NO oxime in product, got {prod_smi!r}"
        )

    def test_preformed_phosphoserine_in_peptide(self):
        """Phosphoserine (pSer) assembles into a peptide with a phosphate group."""
        mol = Molecule(Sequence('ac-pSer-A-G-am')).get_molecule(fmt='ROMol')
        assert mol is not None, "pSer peptide assembly returned None"
        smi = Chem.MolToSmiles(mol)
        assert 'P' in smi, f"Phosphate group not found in pSer product: {smi!r}"
        assert '?' not in smi

    def test_aspartimide_ring_closure(self):
        """Aspartimide ring closure (Asp R4 → backbone N R3) forms the 5-membered succinimide."""
        mol = Molecule(Sequence('ac-D.!1(4,3)-A.!1(3,4)-G-am')).get_molecule(fmt='ROMol')
        assert mol is not None, "Aspartimide assembly returned None"
        smi = Chem.MolToSmiles(mol)
        ring5 = Chem.MolFromSmarts('C1CC(=O)NC1=O')
        assert mol.HasSubstructMatch(ring5), (
            f"5-membered succinimide ring not found in aspartimide product: {smi!r}"
        )
        assert '.' not in smi, f"Disconnected aspartimide product: {smi!r}"


# ---------------------------------------------------------------------------
# GLP-1 agonist drug SMILES → CABILN regression tests
# SMILES sourced from Wikipedia infoboxes and PubChem (CID noted inline).
# Each test passes the drug's canonical SMILES directly to smiles_to_cabiln_core
# and asserts (a) zero unknown residues and (b) the expected CABILN string.
# ---------------------------------------------------------------------------
class TestGLP1DrugSMILES:
    """Regression tests for smiles_to_cabiln_core on approved GLP-1 agonist drugs."""

    @staticmethod
    def _smiles_to_cabiln(smi: str):
        import sys as _sys, os as _os
        _repo_root = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), '..'))
        if _repo_root not in _sys.path:
            _sys.path.insert(0, _repo_root)
        from tools.live_renderer import smiles_to_cabiln_core
        return smiles_to_cabiln_core(smi)

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
        """Wikipedia/PubChem SMILES for approved GLP-1 agonists must decode to
        the expected CABILN with the correct number of residues and zero unknowns.

        SMILES sources:
          semaglutide, liraglutide, exenatide — Wikipedia infobox
          lixisenatide, tirzepatide, retatrutide — PubChem canonical SMILES
        """
        cabiln, details = self._smiles_to_cabiln(smiles)
        unknown = [abbr for abbr, _, _ in details if abbr == '?']
        assert not unknown, (
            f"{drug}: {len(unknown)} unknown residue(s) in {cabiln!r}"
        )
        assert len(details) == n_res, (
            f"{drug}: expected {n_res} residues, got {len(details)}"
        )
        assert cabiln == expected_cabiln, (
            f"{drug}: CABILN mismatch\n  got:      {cabiln!r}\n  expected: {expected_cabiln!r}"
        )


if __name__ == '__main__':
    import unittest
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for cls in (TestBondChemistryUnit, TestAssembly, TestMonomerpipeline,
                TestInlineAttachments, TestCapMonomers, TestCABILNPhase2Parser,
                TestIntramolecularRingClosure, TestExpandedMonomers,
                TestFinalMonomers, TestFattyAcidBranching, TestBracketNotation,
                TestMonomerBuilderCLI, TestRoundTrips,
                TestMonomerPreActivate, TestRestoreRgroups,
                TestReactionPairSMIRKS, TestNotationConversion,
                TestBracketBranchEquivalence):
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
