"""Bundled monomer chemistry, aliases, and reviewed library audit."""

import os
import sys

import pytest
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


from _chemistry_oracles import (
    _romol,
    _smiles,
    _natoms,
)


class TestExpandedMonomers:
    """Smoke-tests for monomers migrated from the legacy 322-entry library."""

    @pytest.mark.parametrize(
        'notation,expected',
        [
            pytest.param('dW-A', 20, id='d_tryptophan_dipeptide'),
            pytest.param('dP-A', 13, id='d_proline_dipeptide'),
            pytest.param('Nle-A', 14, id='norleucine_dipeptide'),
            pytest.param('Hyp-A', 14, id='hydroxyproline_dipeptide'),
            pytest.param('A-Aib-A', 17, id='aib_tripeptide'),
            pytest.param('E_g-A', 15, id='gamma_glutamic_acid_dipeptide'),
            pytest.param('aMePhe-A', 18, id='aMePhe_dipeptide'),
            pytest.param('1Nal-A', 21, id='1_naphthylalanine_dipeptide'),
            pytest.param('2Nal-A', 21, id='2_naphthylalanine_dipeptide'),
            pytest.param('Cha-A', 17, id='cyclohexylalanine_dipeptide'),
            pytest.param('A-dA-A', 16, id='mixed_d_l_tripeptide'),
            pytest.param('Nle-Cha-A', 25, id='nle_cha_tripeptide'),
            pytest.param('Cbz-A-am', 16, id='cbz_cap_assembly'),
            pytest.param('Boc-A-am', 13, id='boc_cap_assembly'),
        ],
    )
    def test_product_atom_count(self, notation, expected):
        assert _natoms(notation) == expected

    @pytest.mark.parametrize(
        'notation,minimum',
        [
            pytest.param('dA-A', 10, id='d_alanine_dipeptide'),
            pytest.param('dF-A', 16, id='d_phenylalanine_dipeptide'),
            pytest.param('Nva-A', 11, id='norvaline_dipeptide'),
            pytest.param('Orn-A', 12, id='ornithine_dipeptide'),
            pytest.param('Hse-A', 11, id='homoserine_dipeptide'),
            pytest.param('Hcy-A', 11, id='homocysteine_dipeptide'),
            pytest.param('Phe_4F-A', 16, id='phe_4F_dipeptide'),
            pytest.param('Phe_4Cl-A', 16, id='phe_4Cl_dipeptide'),
        ],
    )
    def test_product_minimum_atom_count(self, notation, minimum):
        assert _natoms(notation) > minimum


class TestFinalMonomers:
    """Tests for the last 31 monomers completing the 322-entry library."""

    @pytest.mark.parametrize(
        'notation,expected',
        [
            pytest.param('G-Ala_ol', 9, id='ala_ol'),
            pytest.param('G-Gly_ol', 7, id='gly_ol'),
            pytest.param('G-Pro_ol', 11, id='pro_ol'),
            pytest.param('G-Val_ol', 11, id='val_ol'),
            pytest.param('G-Leu_ol', 12, id='leu_ol'),
            pytest.param('G-Phe_ol', 15, id='phe_ol'),
            pytest.param('G-Phg_ol', 14, id='phg_ol'),
            pytest.param('G-Lys_ol', 13, id='lys_ol'),
            pytest.param('G-Thr_ol', 11, id='thr_ol'),
            pytest.param('G-Aib_ol', 10, id='aib_ol'),
            pytest.param('G-D_Phg_ol', 14, id='d_phg_ol'),
            pytest.param('G-D_Pro_ol', 11, id='d_pro_ol'),
            pytest.param('G-D_Thr_ol', 11, id='d_thr_ol'),
            pytest.param('G-Hsl', 11, id='hsl'),
            pytest.param('G-Ala_al', 9, id='ala_al'),
            pytest.param('G-Gly_al', 8, id='gly_al'),
            pytest.param('G-Pro_al', 11, id='pro_al'),
            pytest.param('G-Leu_al', 12, id='leu_al'),
            pytest.param('G-Phe_al', 15, id='phe_al'),
            pytest.param('G-Lys_al', 13, id='lys_al'),
            pytest.param('G-Arg_al', 15, id='arg_al'),
            pytest.param('ac-NMebAla-am', 10, id='nme_beta_ala'),
            pytest.param('ac-NMe2Abz-am', 15, id='nme2abz'),
            pytest.param('ac-NMe4Abz-am', 14, id='nme4abz'),
            pytest.param('NMe23Abz-G-am', 16, id='nme23abz_cap'),
            pytest.param('NMe24Abz-G-am', 16, id='nme24abz_cap'),
            pytest.param('ac-pnA-am', 24, id='pna_adenine'),
            pytest.param('ac-pnC-am', 22, id='pna_cytosine'),
            pytest.param('ac-pnG-am', 25, id='pna_guanine'),
            pytest.param('ac-pnT-am', 23, id='pna_thymine'),
            pytest.param('ac-Pqa-am', 24, id='pqa'),
        ],
    )
    def test_product_atom_count(self, notation, expected):
        assert _natoms(notation) == expected

    def test_pna_tetranucleotide(self):
        """ac-pnA-pnC-pnG-pnT-am: four-base PNA strand."""
        mol = _romol('ac-pnA-pnC-pnG-pnT-am')
        assert mol.GetNumAtoms() > 60


class TestFattyAcidBranching:
    """Sidechain fatty-acid branching via the new amine_primary+backbone_c amide
    reaction."""

    @pytest.mark.parametrize(
        'notation,expected',
        [
            pytest.param('ac-AEEA-am', 14, id='aeea_linear'),
            pytest.param('ac-aMeLeu-am', 13, id='ameLeu_linear'),
            pytest.param('C20FA-E_g-am', 33, id='c20fa_E_g_branch_linear'),
            pytest.param('C20FA-E_g-AEEA-am', 43, id='c20fa_E_g_aeea_linear'),
            pytest.param('ac-K.!1(4,2)-G-am%C20FA-E_g-AEEA-!1', 59,
                id='lys_sidechain_amide_bond'),
        ],
    )
    def test_product_atom_count(self, notation, expected):
        assert _natoms(notation) == expected

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


class TestLibraryRoundTrip:
    """Every record must match the explicit versioned compatibility baseline."""

    def test_library_activation_matches_reviewed_records(self):
        from pyPept.library_quality import audit_library, quality_manifest

        expected = quality_manifest()
        actual = audit_library()
        assert actual["library_sha256"] == expected["library_sha256"], (
            "Bundled definitions changed; review every affected audit record"
        )
        assert actual["records"].keys() == expected["records"].keys()
        changed = [
            symbol for symbol in actual["records"]
            if actual["records"][symbol] != expected["records"][symbol]
        ]
        assert not changed, (
            "Review activation/metadata changes for these exact definitions: "
            + ", ".join(changed)
        )
        assert actual["summary"] == expected["summary"]
