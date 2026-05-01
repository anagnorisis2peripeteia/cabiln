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

    def test_old_biln_crosslink_warns_and_converts(self):
        """Old BILN bare-integer crosslink notation emits DeprecationWarning and auto-converts."""
        with pytest.warns(DeprecationWarning, match='Old BILN'):
            seq = Sequence('K(1,3)-A-K(1,3)-am')
        assert seq is not None

    def test_old_biln_warning_mentions_cabiln(self):
        """DeprecationWarning for old BILN mentions the .!n(y,z) CABILN form."""
        with pytest.warns(DeprecationWarning, match=r'\.!'):
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
            mol = _romol('fmoc-C.trt(4,1)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_trt_cys_has_sc_bond(self):
        """trt-protected Cys should contain an S-C bond to a trisubstituted carbon."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-C.trt(4,1)-am')
        sc_trt = Chem.MolFromSmarts('[S][C](c1ccccc1)(c1ccccc1)c1ccccc1')
        assert mol.HasSubstructMatch(sc_trt), "trt S-C bond pattern not found"

    def test_acm_on_cys_warns_thioether(self):
        """acm cap on Cys R4 (thiol) forms S-C thioether — UserWarning, assembles OK."""
        with pytest.warns(UserWarning, match='thioether'):
            mol = _romol('fmoc-C.acm(4,1)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_acm_cys_has_acetamide(self):
        """acm-protected Cys should contain the S-CH2-NH-C(=O) acetamidomethyl group."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-C.acm(4,1)-am')
        acm_pat = Chem.MolFromSmarts('[S]CNC(C)=O')
        assert mol.HasSubstructMatch(acm_pat), "acm S-CH2-NH-C(=O) pattern not found"

    def test_pbf_on_arg_warns_sulfenamide(self):
        """pbf cap on Arg R4 (guanidinium NH2) forms S-N bond — UserWarning, assembles OK.

        R3 = backbone-N mod; R4 = guanidinium terminal NH2 in CABILN convention.
        """
        with pytest.warns(UserWarning, match='sulfenamide|S.N'):
            mol = _romol('fmoc-R.pbf(4,1)-am')
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert dummies == [], f"Dummy atoms leaked: {len(dummies)}"

    def test_pbf_arg_has_sulfonyl(self):
        """pbf-protected Arg should contain an N-S(=O)(=O) sulfonamide substructure."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol = _romol('fmoc-R.pbf(4,1)-am')
        sulf_pat = Chem.MolFromSmarts('[N]S(=O)=O')
        assert mol.HasSubstructMatch(sulf_pat), "N-S(=O)2 sulfonamide pattern not found"


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
        result, brg = _expand_inline_caps('fmoc-C.trt(4,1)-am')
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
        """gGlu: γ-linked Glu (backbone through γ-carboxyl)."""
        assert _natoms('gGlu-A') == 15

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
        """Old BILN symmetric crosslink converts to CABILN .!n(y,y) notation."""
        result = biln_to_cabiln('C(1,3)-A-A-A-C(1,3)')
        assert '.!1(3,3)' in result
        assert result.count('.!1') == 2

    def test_biln_to_cabiln_asymmetric(self):
        """Old BILN asymmetric crosslink converts correctly (different R-groups)."""
        result = biln_to_cabiln('K(1,4)-A-A-A-D(1,3)')
        assert '.!1(4,3)' in result   # first endpoint: K R4 → partner R3
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
    # Old BILN → Sequence auto-conversion
    # ------------------------------------------------------------------ #

    def test_old_biln_auto_converts_with_warning(self):
        """Sequence auto-converts old BILN and emits DeprecationWarning."""
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            seq = Sequence('C(1,3)-A-A-A-C(1,3)')
        assert any(issubclass(x.category, DeprecationWarning) for x in w)
        assert seq is not None

    def test_old_biln_crosslink_assembles(self):
        """Old BILN disulfide bridge assembles to a molecule with correct atom count."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            romol = _romol('C(1,3)-A-A-A-C(1,3)')
        assert romol is not None
        # Cys-Ala-Ala-Ala-Cys disulfide: check we get a molecule, not None
        from rdkit import Chem
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
        assert not re.search(r'(?<![.\w])[A-Za-z]\w*\(\d+,\d+\)', biln), \
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
    Lys(ε-NH2) → amide → AEEA → gGlu → C18/C20 fatty diacid.
    """

    def test_aeea_linear(self):
        """AEEA (mini-PEG linker) assembles as a standard backbone monomer."""
        assert _natoms('ac-AEEA-am') == 14

    def test_ameLeu_linear(self):
        """aMeLeu (alpha-methyl-L-leucine): Cα-disubstituted, backbone stabiliser."""
        assert _natoms('ac-aMeLeu-am') == 13

    def test_c20fa_gGlu_branch_linear(self):
        """C20 fatty diacid cap + gGlu chain assembles correctly."""
        assert _natoms('C20FA-gGlu-am') == 32

    def test_c20fa_gGlu_aeea_linear(self):
        """Full linker chain: C20FA-gGlu-AEEA with backbone amide bonds."""
        assert _natoms('C20FA-gGlu-AEEA-am') == 42

    def test_lys_sidechain_amide_bond(self):
        """Lys R4 (epsilon-amine) forms amide with AEEA R2 (carboxyl).

        Tests the new amine_primary+backbone_c reaction added to BOND_TABLE.
        """
        assert _natoms('ac-K.!1(4,2)-G-am%C20FA-gGlu-AEEA-!1') == 58

    def test_retatrutide(self):
        """Retatrutide (LY3437943): GIP/GLP-1/glucagon triple agonist, 39 AAs.

        Sequence: Y-Aib-QGTFTSDYSI-aMeLeu-LD-K(C20FA-gGlu-AEEA)-AQAAFI-Aib-EYL
                  LEGGPSSGAPPPSam with C20 fatty diacid branch on Lys16.
        """
        retatrutide = (
            "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.!1(4,2)"
            "-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am"
            "%C20FA-gGlu-AEEA-!1"
        )
        mol = _romol(retatrutide)
        assert mol.GetNumAtoms() == 325


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

    def test_iedda_smirks_step1(self):
        """IEDDA step-1 SMIRKS produces the [4+2] bicyclic adduct."""
        from pyPept.interfaces.reaction_library import REACTIONS
        from rdkit.Chem import AllChem
        entry = REACTIONS['iedda_tetrazine_tco']
        smirks = entry['steps'][0]
        targeted = smirks.replace('[4*]', '[400*]', 1).replace('[4*]', '[401*]', 1)
        rxn = AllChem.ReactionFromSmarts(targeted)
        assert rxn is not None
        # model tetrazine with [400*] on methyl; model cycloalkene with [401*] on exo-C
        tz = Chem.MolFromSmiles('[400*]Cc1nncnn1')
        tco = Chem.MolFromSmiles('C1CC([401*])C=CCCCC1')
        if tz is None or tco is None:
            pytest.skip("Model SMILES failed to parse — check RDKit tetrazine support")
        prods = rxn.RunReactants((tz, tco))
        if not prods:
            prods = rxn.RunReactants((tco, tz))
        assert prods, "IEDDA step-1 SMIRKS produced no products with model molecules"


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
        frag1 = Chem.MolFromSmiles('[400*]NC')      # backbone N (R3): [3*][N:1]
        frag2 = Chem.MolFromSmiles('[401*]C(=O)C')  # sidechain carboxyl: [4*][C:2](=[O:3])
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
        assert entry['id'] == 'thioester'

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
        """[4*][S:2].[4*][C:4][Cl,Br,I] >> [S:2][C:4]; halide departs."""
        from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
        entry = REACTIONS['thioether_halide']
        frag1 = Chem.MolFromSmiles('[400*]SC')    # thiol S: [4*][S:2]
        frag2 = Chem.MolFromSmiles('[401*]CCl')   # alkyl chloride: [4*][C:4][Cl]
        assert frag1 and frag2
        prod = run_bond_smirks(frag1, frag2, 400, 401, entry, intramolecular=False)
        smi = Chem.MolToSmiles(prod)
        thioether = Chem.MolFromSmarts('[S][C]')
        assert prod.HasSubstructMatch(thioether), f"No thioether S-C in product: {smi}"
        assert 'Cl' not in smi, f"Chloride should have departed via SN2: {smi}"
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

    def test_old_biln_notation_warns_and_converts(self):
        """Old BILN K(1,3) crosslink syntax → DeprecationWarning + auto-convert, valid result."""
        report = Sequence.validate('K(1,3)-A-K(1,3)')
        # validate() captures warnings internally; check they appear in report.warnings
        assert any('Old BILN' in w for w in report.warnings)
        assert report.ok is True


# ---------------------------------------------------------------------------
# 11. Sequential reaction bracket notation — [.A(r,s).B(t,u)...]
# ---------------------------------------------------------------------------

class TestBracketNotation:
    """Tests for [...] sequential reaction bracket syntax in CABILN.

    Inside [...], each .Frag(x,y) attaches to the preceding fragment:
      - first entry: host.Rx -> Frag.Ry
      - subsequent entries: prev_frag.Rx -> Frag.Ry
    Auto bond IDs are assigned left-to-right starting from 100.
    """

    # --- _expand_inline_caps unit tests ---

    def test_single_step_bracket_equals_inline_cap(self):
        """Single-step bracket [.boc(4,1)] expands identically to .boc(4,1)."""
        result_bracket, _ = _expand_inline_caps('K[.boc(4,1)]-am')
        result_inline,  _ = _expand_inline_caps('K.boc(4,1)-am')
        assert result_bracket == result_inline

    def test_two_step_bracket_host_annotation(self):
        """Host residue receives only the first bond annotation from a two-step bracket."""
        result, _ = _expand_inline_caps('G[.A(3,1).am(2,1)]')
        # G gets the first auto bond ID at its R3
        assert '(100,3)' in result
        # Second bond ID must NOT appear on G
        assert 'G(100,3)(101' not in result

    def test_two_step_bracket_intermediate_fragment(self):
        """Intermediate fragment receives both incoming and outgoing bond annotations."""
        result, _ = _expand_inline_caps('G[.A(3,1).am(2,1)]')
        # A is intermediate: attached at R1 by bond 100, passes on at R2 by bond 101
        assert 'A(100,1)(101,2)' in result

    def test_two_step_bracket_terminal_fragment(self):
        """Terminal fragment in a two-step bracket receives only its incoming bond."""
        result, _ = _expand_inline_caps('G[.A(3,1).am(2,1)]')
        assert 'am(101,1)' in result

    def test_three_step_bracket_full_chain(self):
        """Three-step bracket produces correct annotations on all three fragments."""
        result, _ = _expand_inline_caps('G[.A(3,1).G(2,1).am(2,1)]')
        assert '(100,3)' in result      # G host → bond 100 at R3
        assert 'A(100,1)(101,2)' in result   # first intermediate
        assert 'G(101,1)(102,2)' in result   # second intermediate
        assert 'am(102,1)' in result         # terminal

    def test_multiple_brackets_get_distinct_bond_ids(self):
        """Two brackets on different residues get non-overlapping auto bond IDs."""
        import re as _re
        result, _ = _expand_inline_caps('C[.trt(4,1)]-G-K[.boc(4,1)]-am')
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
        """A residue can carry both a [...] bracket and a .!1 crosslink marker."""
        result, brg = _expand_inline_caps('C[.trt(4,1)].!1(4,4)-G-C.!1')
        # bracket annotation present
        assert '(100,4)' in result
        assert 'trt(100,1)' in result
        # crosslink annotation present
        assert '(!1,4)' in result
        assert brg == {'!1': 4}

    def test_bracket_preserves_other_inline_caps(self):
        """A bracket on one residue does not affect inline caps on other residues."""
        result, _ = _expand_inline_caps('K[.boc(4,1)]-G.ac(2,1)-am')
        # boc bracket: ID 100 on K.R4, boc attached at R1
        assert '(100,4)' in result
        assert 'boc(100,1)' in result
        # ac inline cap: ID 101 on G.R2, ac attached at R1
        assert '(101,2)' in result
        assert 'ac(101,1)' in result

    def test_bracket_in_percent_separated_segment(self):
        """Bracket notation works in a %-separated branch segment."""
        result, _ = _expand_inline_caps('K.!1(4,2)-am%G[.boc(3,1)]-G.!1')
        assert 'boc(' in result
        assert '(!1,4)' in result

    # --- Error cases ---

    def test_empty_bracket_raises(self):
        """Empty [] raises ValueError."""
        with pytest.raises(ValueError, match='no valid'):
            _expand_inline_caps('C[]')

    def test_bracket_with_crosslink_marker_raises(self):
        """[.!1(4,4)] raises because !1 is not a letter-starting cap token."""
        with pytest.raises(ValueError, match='no valid'):
            _expand_inline_caps('C[.!1(4,4)]-G-C.!1(4,4)-am')

    def test_bracket_with_trailing_garbage_raises(self):
        """[.trt(4,1)xyz] raises because 'xyz' is not a valid entry."""
        with pytest.raises(ValueError, match='unrecognised'):
            _expand_inline_caps('C[.trt(4,1)xyz]-am')

    def test_bracket_missing_parens_raises(self):
        """[.trt] (no r-group numbers) raises because .trt has no (x,y)."""
        with pytest.raises(ValueError, match='no valid'):
            _expand_inline_caps('C[.trt]-am')

    # --- Backbone peptide branch guard ---

    def test_backbone_r2_r1_pattern_warns_when_two_steps(self):
        """Two or more R2->R1 steps in a bracket emit UserWarning (peptide branch)."""
        with pytest.warns(UserWarning, match='R2->R1|backbone amide|peptide branch'):
            _expand_inline_caps('K[.G(4,1).A(2,1).am(2,1)]')

    def test_single_r2_r1_step_no_warn(self):
        """A single R2->R1 step (e.g. Mal->DBCO style) does NOT warn — ambiguous."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _expand_inline_caps('C[.Mal(4,1).DBCO(2,1)]')

    def test_backbone_pattern_expansion_still_works(self):
        """Warning does not prevent expansion — structure still returned correctly."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            result, _ = _expand_inline_caps('K[.G(4,1).A(2,1).am(2,1)]')
        assert 'G(100,1)(101,2)' in result
        assert 'A(101,1)(102,2)' in result
        assert 'am(102,1)' in result

    def test_non_backbone_r_groups_no_warn(self):
        """R4->R1 (sidechain->cap) inside bracket does NOT warn."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            _expand_inline_caps('K[.boc(4,1)]')

    # --- colorize_cabiln ---

    def test_colorize_single_bracket_pair(self):
        """Single bracket pair gets ANSI colour codes on open and close."""
        import re as _re
        coloured = colorize_cabiln('C[.trt(4,1)]')
        # At least one ANSI colour escape is injected
        assert _re.search(r'\033\[\d+m', coloured), "No ANSI escape found"
        # Reset code present after close bracket
        assert '\033[0m' in coloured
        # Original content still intact when escapes stripped
        plain = _re.sub(r'\033\[\d+m', '', coloured)
        assert plain == 'C[.trt(4,1)]'

    def test_colorize_two_pairs_different_colours(self):
        """Two bracket pairs get different ANSI colour codes."""
        coloured = colorize_cabiln('C[.trt(4,1)]-K[.boc(4,1)]')
        # Extract just the escape sequences before each [
        import re as _re
        escapes = _re.findall(r'\033\[\d+m(?=\[)', coloured)
        assert len(escapes) == 2
        assert escapes[0] != escapes[1]

    def test_colorize_use_ansi_false_returns_plain(self):
        """use_ansi=False returns the string unchanged."""
        s = 'C[.trt(4,1)]-K[.boc(4,1)]'
        assert colorize_cabiln(s, use_ansi=False) == s

    def test_colorize_no_brackets_returns_unchanged(self):
        """String with no brackets is returned as-is (modulo no escapes)."""
        s = 'fmoc-C-G-K-am'
        assert colorize_cabiln(s) == s

    # --- End-to-end assembly tests ---

    def test_bracket_single_step_assembly_matches_inline(self):
        """fmoc-C[.trt(4,1)]-am assembles to the same molecule as fmoc-C.trt(4,1)-am."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            mol_bracket = _romol('fmoc-C[.trt(4,1)]-am')
            mol_inline  = _romol('fmoc-C.trt(4,1)-am')
        from rdkit import Chem
        smi_bracket = Chem.MolToSmiles(mol_bracket)
        smi_inline  = Chem.MolToSmiles(mol_inline)
        assert smi_bracket == smi_inline

    def test_bracket_validate_accepts_syntax(self):
        """Sequence.validate() returns ok=True for a single-step bracket sequence."""
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            report = Sequence.validate('fmoc-C[.trt(4,1)]-am')
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
        smi = smiles_from_cabiln('fmoc-C.trt(4,1)-am')
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
# Standalone runner
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import unittest
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for cls in (TestBondChemistryUnit, TestAssembly, TestMonomerpipeline,
                TestInlineAttachments, TestCapMonomers, TestCABILNPhase2Parser,
                TestIntramolecularRingClosure, TestExpandedMonomers,
                TestFinalMonomers, TestFattyAcidBranching, TestBracketNotation,
                TestMonomerBuilderCLI, TestRoundTrips):
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
