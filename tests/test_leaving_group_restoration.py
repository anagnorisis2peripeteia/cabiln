"""Leaving-group restoration and products from activated monomers."""

import os
import sys

import pytest
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyPept.molecule import Molecule
from pyPept.sequence import Sequence


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

    @pytest.mark.filterwarnings("ignore")
    @pytest.mark.parametrize(
        'notation,expected',
        [
            pytest.param('ac-DGlu-am', {'acid': 1, 'dummy': 0},
                id='d_hasp_carboxyl_restores'),
            pytest.param('ac-hGlu-am', {'acid': 1, 'dummy': 0},
                id='hglu_carboxyl_restores'),
            pytest.param('ac-D_hGlu-am', {'acid': 1, 'dummy': 0},
                id='d_hglu_carboxyl_restores'),
            pytest.param('ac-b3hAsp-am', {'acid': 1, 'dummy': 0},
                id='b3hasp_carboxyl_restores'),
            pytest.param('ac-D_b3hAsp-am', {'acid': 1, 'dummy': 0},
                id='d_b3hasp_carboxyl_restores'),
            pytest.param('ac-b3hGlu-am', {'acid': 1, 'dummy': 0},
                id='b3hglu_carboxyl_restores'),
            pytest.param('ac-D_b3hGlu-am', {'acid': 1, 'dummy': 0},
                id='d_b3hglu_carboxyl_restores'),
            pytest.param('ac-aMeAsp-am', {'acid': 1, 'dummy': 0},
                id='ameasp_carboxyl_restores'),
            pytest.param('ac-aMeGlu-am', {'acid': 1, 'dummy': 0},
                id='ameglu_carboxyl_restores'),
            pytest.param('ac-D_meD-am', {'acid': 1, 'dummy': 0},
                id='d_med_slot3_carboxyl_restores'),
            pytest.param('ac-D_meE-am', {'acid': 1, 'dummy': 0},
                id='d_mee_slot3_carboxyl_restores'),
            pytest.param('ac-D-am', {'acid': 1, 'dummy': 0},
                id='asp_c_attach_carboxyl_restores'),
            pytest.param('ac-E-am', {'acid': 1, 'dummy': 0},
                id='glu_c_attach_carboxyl_restores'),
            pytest.param('ac-G-am', {'acid': 0, 'nh2': 0, 'dummy': 0},
                id='both_caps_suppress_both_termini'),
            pytest.param('ac-C-am', {'thiol': 1, 'acid': 0, 'dummy': 0},
                id='free_thiol_restores_on_cysteine'),
            pytest.param('ac-K-am', {'nh2': 1, 'acid': 0, 'dummy': 0},
                id='free_amine_restores_on_lysine'),
            pytest.param('ac-D-D-D-am', {'acid': 3, 'dummy': 0},
                id='triple_asp_three_carboxyls_restored'),
            pytest.param('ac-D-D-D-D-D-am', {'acid': 5, 'dummy': 0},
                id='five_asp_five_carboxyls_restored'),
            pytest.param('ac-E-hGlu-am', {'acid': 2, 'dummy': 0},
                id='mixed_hasp_hglu_two_carboxyls'),
            pytest.param('ac-D-E-am', {'acid': 2, 'dummy': 0},
                id='asp_glu_two_carboxyls_restored'),
            pytest.param('ac-C-D-am', {'thiol': 1, 'acid': 1, 'dummy': 0},
                id='cys_and_asp_independent_restore'),
            pytest.param('ac-C.!1(4,4)-D-C.!1-am',
                {'thiol': 0, 'acid': 1, 'dummy': 0},
                id='disulfide_eliminates_thiols_but_asp_cooh_survives'),
            pytest.param('ac-K.!1(4,4)-G-D.!1-am', {'acid': 0, 'dummy': 0},
                id='lactam_crosslink_consumes_carboxyl_not_restored'),
            pytest.param('A-G-A', {'dummy': 0}, id='no_dummies_simple_tripeptide'),
            pytest.param('ac-D-K-C-E-am',
                {'acid': 2, 'thiol': 1, 'nh2': 1, 'dummy': 0},
                id='no_dummies_complex_mixed_sequence'),
            pytest.param('ac-C.!1(4,4)-A-C.!1-am', {'dummy': 0},
                id='no_dummies_after_disulfide_crosslink'),
            pytest.param('ac-S-am', {'hydroxyl': 1, 'acid': 0, 'dummy': 0},
                id='ser_hydroxyl_restores'),
            pytest.param('ac-T-am', {'hydroxyl': 1, 'acid': 0, 'dummy': 0},
                id='thr_hydroxyl_restores'),
            pytest.param('ac-Y-am', {'phenol': 1, 'dummy': 0},
                id='tyr_phenolic_oh_restores'),
            pytest.param('ac-W-am', {'ar_nh': 1, 'dummy': 0},
                id='trp_indole_nh_restores'),
            pytest.param('ac-H-am', {'ar_nh': 1, 'dummy': 0},
                id='his_imidazole_nh_restores'),
            pytest.param('ac-M-am', {'thiol': 0, 'dummy': 0},
                id='met_thioether_not_thiol'),
            pytest.param('G', {'nh2': 1, 'acid': 1, 'dummy': 0},
                id='single_gly_free_termini'),
            pytest.param('D', {'nh2': 1, 'acid': 2, 'dummy': 0},
                id='single_asp_free_termini'),
            pytest.param('E', {'nh2': 1, 'acid': 2, 'dummy': 0},
                id='single_glu_free_termini'),
            pytest.param('C', {'nh2': 1, 'acid': 1, 'thiol': 1, 'dummy': 0},
                id='single_cys_free_termini'),
            pytest.param('A-G-D', {'nh2': 1, 'acid': 2, 'dummy': 0},
                id='free_n_terminus_tripeptide'),
            pytest.param('ac-D', {'nh2': 0, 'acid': 2, 'dummy': 0},
                id='ac_only_suppresses_nh2'),
            pytest.param('D-am', {'nh2': 1, 'acid': 1, 'dummy': 0},
                id='am_only_suppresses_c_term_acid'),
            pytest.param('ac-D-D-D-D-D-D-D-D-D-D-am', {'acid': 10, 'dummy': 0},
                id='ten_asp_ten_carboxyls'),
            pytest.param('ac-C-C-C-C-C-am', {'thiol': 5, 'dummy': 0},
                id='five_cys_five_thiols'),
            pytest.param('ac-C-D-K-am',
                {'thiol': 1, 'acid': 1, 'nh2': 1, 'dummy': 0},
                id='mixed_acid_thiol_amine'),
            pytest.param('ac-S-D-C-am',
                {'hydroxyl': 1, 'acid': 1, 'thiol': 1, 'dummy': 0},
                id='mixed_hydroxyl_acid_thiol'),
            pytest.param('ac-W-D-am', {'ar_nh': 1, 'acid': 1, 'dummy': 0},
                id='mixed_aromatic_nh_acid'),
            pytest.param('ac-A-A-A-am',
                {'acid': 0, 'thiol': 0, 'nh2': 0, 'dummy': 0},
                id='all_ala_no_sidechain_restore'),
            pytest.param('ac-G-G-G-am',
                {'acid': 0, 'thiol': 0, 'nh2': 0, 'dummy': 0},
                id='all_gly_no_sidechain_restore'),
            pytest.param('ac-b3hAsp-b3hGlu-am', {'acid': 2, 'dummy': 0},
                id='beta_homo_asp_glu_mix'),
            pytest.param('ac-b3hAsp-D-am', {'acid': 2, 'dummy': 0},
                id='beta_and_standard_asp_mix'),
            pytest.param('ac-aMeAsp-E-am', {'acid': 2, 'dummy': 0},
                id='ameasp_and_standard_glu_mix'),
            pytest.param('ac-K.!1(4,4)-D.!1-D-am', {'acid': 1, 'dummy': 0},
                id='lactam_consumes_one_asp_other_survives'),
            pytest.param('ac-C.!1(4,4)-D-K-C.!1-am',
                {'thiol': 0, 'acid': 1, 'nh2': 1, 'dummy': 0},
                id='disulfide_crosslink_acid_and_amine_survive'),
        ],
    )
    def test_restored_group_counts(self, notation, expected):
        actual = self._counts(notation)
        for group, count in expected.items():
            assert actual[group] == count, (notation, group, actual)

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
