"""Bond validation, molecule assembly, and assembly error contracts."""

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
    _check_bond_chemistry,
)

from _chemistry_oracles import (
    _romol,
    _smiles,
    _natoms,
    _assert_same_monomer_partition,
)


def _mol(smi):
    m = Chem.MolFromSmiles(smi)
    assert m is not None, f"Bad SMILES: {smi}"
    return m


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

    def test_sanitize_error_wrapped_as_valueerror(self, monkeypatch):
        """A sanitization failure carries attachment context and its RDKit cause."""
        from pyPept.assembly_cache import clear_assembly_cache

        clear_assembly_cache()
        sequence = Sequence('A')

        def invalid_valence(_molecule):
            raise ValueError('invalid valence control')

        monkeypatch.setattr(Chem, 'SanitizeMol', invalid_valence)
        with pytest.raises(ValueError, match='R-group assignments') as error:
            Molecule(sequence, depiction=None)
        assert str(error.value.__cause__) == 'invalid valence control'

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


class TestCapMonomers:
    """Supported protecting groups retain their chemistry without false alerts."""

    @pytest.mark.parametrize('source,expected_bond', [
        ('fmoc-C.trt(4,2)-am', '[S][C](c1ccccc1)(c1ccccc1)c1ccccc1'),
        ('fmoc-C.acm(4,2)-am', '[S]CNC(C)=O'),
        ('fmoc-R.pbf(5,2)-am', '[N]S(=O)=O'),
        ('G.OBn_(1,1)-am', 'c1ccccc1CONCC(N)=O'),
        ('G.OMe_(1,1)-am', 'CONCC(N)=O'),
    ])
    def test_supported_cap_assembles_without_a_chemistry_warning(self, source, expected_bond):
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            mol = _romol(source)
        assert not any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        assert mol.HasSubstructMatch(Chem.MolFromSmarts(expected_bond))


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


class TestEdgeCasesErrorHandling:
    """Tests for boundary conditions and invalid-input behaviour."""

    @staticmethod
    def _setup():
        from pyPept.smiles import smiles_to_cabiln_core
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


class TestAdditionalChemistryEdgeCases:
    """Assembly and SMIRKS correctness for less-common reaction types."""

    def test_thiol_maleimide_adduct_assembly(self):
        """Thiol-maleimide ligation (Cys + 5FM) forms the expected thioether-succinimide."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mol = Molecule(Sequence('ac-C.!1(4,6)-am%5FM.!1(6,4)')).get_molecule(fmt='ROMol')
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

    @pytest.mark.parametrize('cabiln,description', [
        ('ac-D_Lac-A-am',     'D-lactic acid N-terminal depsi'),
        ('ac-L_Lac-A-am',     'L-lactic acid N-terminal depsi'),
        ('ac-GlyAc-A-am',     'glycolic acid N-terminal depsi'),
        ('ac-D_Lac-G-am',     'D-lactic acid + Gly'),
        ('ac-D_Lac-A-G-am',   'D-lactic acid + 2 residues'),
        ('ac-A-D_Lac-A-am',   'D-lactic acid mid-chain'),
        ('ac-D_Lac-D_Lac-am', 'two consecutive depsi residues'),
    ])
    def test_depsipeptide_backbone_ester(self, cabiln, description):
        """Backbone ester bond (depsipeptide): D_Lac/L_Lac/GlyAc residues assemble and round-trip."""
        from pyPept.smiles import smiles_to_cabiln_core
        ester_pat = Chem.MolFromSmarts('[C:1](=O)[O:2][C:3]')
        mol = Molecule(Sequence(cabiln)).get_molecule(fmt='ROMol')
        assert mol is not None, f"Assembly returned None for {cabiln!r}"
        smi = Chem.MolToSmiles(mol)
        assert '*' not in smi, f"Unreacted dummy in {cabiln!r}: {smi!r}"
        assert '.' not in smi, f"Disconnected product for {cabiln!r}: {smi!r}"
        assert mol.HasSubstructMatch(ester_pat), f"No ester bond in {cabiln!r}: {smi!r}"
        result, _ = smiles_to_cabiln_core(smi)
        _assert_same_monomer_partition(cabiln, result)

    def test_diselenide_assembly(self):
        """Selenocysteine (Sec) diselenoide crosslink assembles with Se–Se bond."""
        mol = Molecule(Sequence('ac-Sec.!1(4,4)-A-A-Sec.!1(4,4)-am')).get_molecule(fmt='ROMol')
        assert mol is not None, "Diselenide assembly returned None"
        smi = Chem.MolToSmiles(mol)
        sese = Chem.MolFromSmarts('[Se][Se]')
        assert mol.HasSubstructMatch(sese), f"Se–Se bond not found in product: {smi!r}"
        assert '*' not in smi, f"Unreacted dummy in diselenide product: {smi!r}"
        assert '.' not in smi, f"Disconnected diselenide product: {smi!r}"

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
