"""Raw input normalization, monomer pre-activation, and CSV ingestion."""

import csv
import os
import pathlib
import sys
import tempfile

import pytest
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


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

    @pytest.mark.parametrize(
        'token,expected',
        [('T', 'C[C@@H](O)[C@H](N)C(=O)O'),
         ('DThr', 'C[C@H](O)[C@@H](N)C(=O)O')],
    )
    def test_threonine_tokens_preserve_both_stereocenters(self, token, expected):
        from pyPept.interfaces.monomer_pipeline import normalize_input

        smiles, is_template, error = normalize_input(token)
        assert error is None and not is_template
        assert Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == Chem.MolToSmiles(
            Chem.MolFromSmiles(expected))

    def test_normalize_unknown_returns_error(self):
        """Unrecognisable input returns a non-None error string."""
        from pyPept.interfaces.monomer_pipeline import normalize_input
        _, _, err = normalize_input('ZZZZNOTREAL')
        assert err is not None


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
        # Preserve the public legacy tuple and index interfaces.
        assert tuple(r) == (r.chuckles, r.leaving, r.chem_types, None)
        assert [r[index] for index in range(4)] == list(r)
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

    @pytest.mark.parametrize(
        'smiles,expect_slots,expect_lg,expect_ct',
        [
            pytest.param('N[C@H](CC(=O)O)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6},
                {4: '[OH]'}, {4: 'carboxyl'}, id='d_asp_carboxyl_on_carbon'),
            pytest.param('N[C@H](CCC(=O)O)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6},
                {4: '[OH]'}, {4: 'carboxyl'}, id='d_glu_carboxyl_on_carbon'),
            pytest.param('N[C@@H](CS)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 16}, {4: '[H]'},
                {4: 'thiol', 1: 'backbone_n', 2: 'backbone_c'},
                id='cys_thiol_on_sulfur'),
            pytest.param('N[C@H](CS)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 16}, {4: '[H]'},
                {4: 'thiol'}, id='d_cys_thiol_on_sulfur'),
            pytest.param('N[C@@H](C[SeH])C(=O)O', {1: 7, 2: 6, 3: 7, 4: 34},
                {4: '[H]'}, {4: 'selenol'}, id='sec_selenol_on_selenium'),
            pytest.param('N[C@@H](CO)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8}, {4: '[H]'},
                {4: 'hydroxyl', 1: 'backbone_n', 2: 'backbone_c'},
                id='ser_hydroxyl_on_oxygen'),
            pytest.param('N[C@@H]([C@@H](O)C)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8},
                {4: '[H]'}, {4: 'hydroxyl'}, id='thr_hydroxyl_on_oxygen'),
            pytest.param('N[C@@H](CCCCN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                None, {4: 'amine_primary', 5: 'amine_secondary'},
                id='lys_amine_primary_on_nitrogen'),
            pytest.param('N[C@@H](CCCN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7}, None,
                {4: 'amine_primary', 5: 'amine_secondary'},
                id='orn_amine_primary_on_nitrogen'),
            pytest.param('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
                {1: 7, 2: 6, 3: 7, 4: 7}, None,
                {4: 'aromatic_nh', 1: 'backbone_n', 2: 'backbone_c'},
                id='trp_aromatic_nh'),
            pytest.param('N[C@@H](Cc1cnc[nH]1)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7},
                None, {4: 'aromatic_nh', 1: 'backbone_n', 2: 'backbone_c'},
                id='his_aromatic_nh'),
            pytest.param('N[C@@H](Cc1ccc(O)cc1)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8},
                None, {4: 'hydroxyl_phenolic'}, id='tyr_phenolic_oh'),
            pytest.param('N[C@@H](CC(=O)N)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                None, {4: 'amide_nh', 5: 'amide_nh'}, id='asn_amide_nh2'),
            pytest.param('N[C@@H](CCC(=O)N)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                None, {4: 'amide_nh', 5: 'amide_nh'}, id='gln_amide_nh2'),
            pytest.param('N[C@@H](C)C(=O)O', {1: 7, 2: 6, 3: 7},
                {1: '[H]', 2: '[OH]', 3: '[H]'},
                {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'},
                id='ala_no_sidechain'),
            pytest.param('NCC(=O)O', {1: 7, 2: 6, 3: 7}, None,
                {1: 'backbone_n', 2: 'backbone_c'}, id='gly_no_sidechain'),
            pytest.param('N[C@@H](C(C)C)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='val_no_sidechain'),
            pytest.param('N[C@@H](CC(C)C)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='leu_no_sidechain'),
            pytest.param('N[C@@H]([C@@H](C)CC)C(=O)O', {1: 7, 2: 6, 3: 7}, None,
                None, id='ile_no_sidechain'),
            pytest.param('N[C@@H](Cc1ccccc1)C(=O)O', {1: 7, 2: 6, 3: 7}, None,
                {1: 'backbone_n', 2: 'backbone_c'}, id='phe_no_sidechain'),
            pytest.param('N[C@@H](CCSC)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='met_thioether_not_thiol'),
            pytest.param('OC(=O)[C@@H]1CCCN1', {1: 7, 2: 6}, None,
                {1: 'backbone_n', 2: 'backbone_c'}, id='pro_no_backbone_n_mod'),
            pytest.param('CC(N)(C)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='aib_alpha_methyl_no_sidechain'),
            pytest.param('N[C@H](C)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='d_ala_same_as_l_ala'),
            pytest.param('N[C@H](CO)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8}, None,
                {4: 'hydroxyl'}, id='d_ser_same_as_l_ser'),
            pytest.param('N[C@H](CCCCN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7}, None,
                {4: 'amine_primary', 5: 'amine_secondary'}, id='d_lys_same_as_l_lys'),
            pytest.param('N[C@@H](CCS)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 16}, {4: '[H]'},
                {4: 'thiol'}, id='homocysteine_thiol'),
            pytest.param('N[C@H](CCS)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 16}, None,
                {4: 'thiol'}, id='d_homocysteine_thiol'),
            pytest.param('N[C@@H](CCCS)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 16}, None,
                {4: 'thiol'}, id='mercaptonorvaline_thiol'),
            pytest.param('N[C@H](C[SeH])C(=O)O', {1: 7, 2: 6, 3: 7, 4: 34},
                {4: '[H]'}, {4: 'selenol'}, id='d_sec_selenol'),
            pytest.param('N[C@@H](CC[SeH])C(=O)O', {1: 7, 2: 6, 3: 7, 4: 34}, None,
                {4: 'selenol'}, id='homoselenocysteine_selenol'),
            pytest.param('N[C@H](CC[SeH])C(=O)O', {1: 7, 2: 6, 3: 7, 4: 34}, None,
                {4: 'selenol'}, id='d_homoselenocysteine_selenol'),
            pytest.param('N[C@@H](CCC[SeH])C(=O)O', {1: 7, 2: 6, 3: 7, 4: 34}, None,
                {4: 'selenol'}, id='long_chain_selenol'),
            pytest.param('N[C@@H](CCl)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6}, {4: '[Cl]'},
                {4: 'alkyl_halide_c'}, id='chloroalanine_alkyl_halide'),
            pytest.param('N[C@@H](CBr)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6}, None,
                {4: 'alkyl_halide_c'}, id='bromoalanine_alkyl_halide'),
            pytest.param('N[C@@H](CI)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6}, None,
                {4: 'alkyl_halide_c'}, id='iodoalanine_alkyl_halide'),
            pytest.param('N[C@H](CCl)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6}, None,
                {4: 'alkyl_halide_c'}, id='d_chloroalanine_alkyl_halide'),
            pytest.param('N[C@@H](CCCl)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 6}, None,
                {4: 'alkyl_halide_c'}, id='chlorohomoalanine_alkyl_halide'),
            pytest.param('N[C@@H](CCC(=O)NN)C(=O)O', None, None, {4: 'hydrazide'},
                id='glu_hydrazide'),
            pytest.param('N[C@@H](CC(=O)NN)C(=O)O', None, None, {4: 'hydrazide'},
                id='asp_hydrazide'),
            pytest.param('N[C@H](CCC(=O)NN)C(=O)O', None, None, {4: 'hydrazide'},
                id='d_glu_hydrazide'),
            pytest.param('N[C@H](CC(=O)NN)C(=O)O', None, None, {4: 'hydrazide'},
                id='d_asp_hydrazide'),
            pytest.param('N[C@@H](CCCC(=O)NN)C(=O)O', None, None, {4: 'hydrazide'},
                id='long_chain_hydrazide'),
            pytest.param('N[C@@H](CCN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7},
                {5: '[H]'}, {4: 'amine_primary', 5: 'amine_secondary'},
                id='dab_amine_primary'),
            pytest.param('N[C@@H](CN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7}, None,
                {4: 'amine_primary', 5: 'amine_secondary'}, id='dap_amine_primary'),
            pytest.param('N[C@H](CCCN)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7, 5: 7}, None,
                {4: 'amine_primary'}, id='d_orn_amine_primary'),
            pytest.param('N[C@H](Cc1ccc(O)cc1)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8},
                None, {4: 'hydroxyl_phenolic'}, id='d_tyr_phenolic_oh'),
            pytest.param('N[C@@H](Cc1cc(F)c(O)cc1)C(=O)O', None, None,
                {4: 'hydroxyl_phenolic'}, id='3f_tyr_phenolic_oh'),
            pytest.param('N[C@@H](Cc1cc(Cl)c(O)cc1)C(=O)O', None, None,
                {4: 'hydroxyl_phenolic'}, id='3cl_tyr_phenolic_oh'),
            pytest.param('N[C@@H](CCO)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8}, None,
                {4: 'hydroxyl'}, id='homoserine_hydroxyl'),
            pytest.param('N[C@H]([C@H](O)C)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 8}, None,
                {4: 'hydroxyl'}, id='d_thr_hydroxyl'),
            pytest.param('N[C@H](Cc1c[nH]c2ccccc12)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7},
                None, {4: 'aromatic_nh'}, id='d_trp_aromatic_nh'),
            pytest.param('N[C@H](Cc1cnc[nH]1)C(=O)O', {1: 7, 2: 6, 3: 7, 4: 7}, None,
                {4: 'aromatic_nh'}, id='d_his_aromatic_nh'),
            pytest.param('N[C@@H](Cc1c[nH]c2cc(F)ccc12)C(=O)O',
                {1: 7, 2: 6, 3: 7, 4: 7}, None, {4: 'aromatic_nh'},
                id='5f_trp_aromatic_nh'),
            pytest.param('N[C@@H](CCC#C)C(=O)O', None, None, {4: 'alkyne_c'},
                id='homopropargylglycine_alkyne'),
            pytest.param('N[C@@H](CCCC#C)C(=O)O', None, None, {4: 'alkyne_c'},
                id='3c_spacer_alkyne'),
            pytest.param('N[C@H](CC#C)C(=O)O', None, None, {4: 'alkyne_c'},
                id='d_propargylglycine_alkyne'),
            pytest.param('N[C@@H](CN=[N+]=[N-])C(=O)O', None, None,
                {4: 'azide_alpha_c'}, id='azidoalanine_azide'),
            pytest.param('N[C@@H](CCN=[N+]=[N-])C(=O)O', None, None,
                {4: 'azide_alpha_c'}, id='azidohomoalanine_azide'),
            pytest.param('N[C@@H](CCCN=[N+]=[N-])C(=O)O', None, None,
                {4: 'azide_alpha_c'}, id='azidoornithine_azide'),
            pytest.param('N[C@@H](CCCCN=[N+]=[N-])C(=O)O', None, None,
                {4: 'azide_alpha_c'}, id='azidolysine_azide'),
            pytest.param('N[C@H](CN=[N+]=[N-])C(=O)O', None, None,
                {4: 'azide_alpha_c'}, id='d_azidoalanine_azide'),
            pytest.param('N[C@@H](CC=C)C(=O)O', None, None, {4: 'terminal_alkene'},
                id='pentenoic_terminal_alkene'),
            pytest.param('N[C@@H](CCC=C)C(=O)O', None, None, {4: 'terminal_alkene'},
                id='hexenoic_terminal_alkene'),
            pytest.param('N[C@H](CCCC=C)C(=O)O', None, None, {4: 'terminal_alkene'},
                id='d_octenoic_terminal_alkene'),
            pytest.param('N[C@@H](CCc1nncnn1)C(=O)O', None, None, {4: 'tetrazine_c'},
                id='tetrazine_2c_spacer'),
            pytest.param('N[C@@H](CCCCCc1nncnn1)C(=O)O', None, None,
                {4: 'tetrazine_c'}, id='tetrazine_5c_spacer'),
            pytest.param('N[C@H](Cc1nncnn1)C(=O)O', None, None, {4: 'tetrazine_c'},
                id='d_tetrazine'),
            pytest.param('N[C@H](CC1=CCCCCCC1)C(=O)O', None, None, {4: 'tco_c'},
                id='d_tco'),
            pytest.param('N[C@@H](CCC1CCCCC#CC1)C(=O)O', None, None,
                {4: 'cyclooctyne_c'}, id='cyclooctyne_2c_linker'),
            pytest.param('N[C@@H](CCCCC1CCCCC#CC1)C(=O)O', None, None,
                {4: 'cyclooctyne_c'}, id='cyclooctyne_4c_linker'),
            pytest.param('N[C@H](CC1CCCCC#CC1)C(=O)O', None, None,
                {4: 'cyclooctyne_c'}, id='d_cyclooctyne'),
            pytest.param('N[C@H](Cc1ccc(C=O)cc1)C(=O)O', None, None, {4: 'aldehyde'},
                id='d_formylphe_aldehyde'),
            pytest.param('N[C@@H](CCC=O)C(=O)O', None, None, {4: 'aldehyde'},
                id='aliphatic_aldehyde'),
            pytest.param('N[C@@H](CCCC=O)C(=O)O', None, None, {4: 'aldehyde'},
                id='long_aliphatic_aldehyde'),
            pytest.param('N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O', None, None,
                {4: 'nhs_ester'}, id='nhs_glu_ester'),
            pytest.param('N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O', None, None,
                {4: 'nhs_ester'}, id='nhs_nva_ester'),
            pytest.param('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O', None, None,
                {4: 'nhs_ester'}, id='d_nhs_asp_ester'),
            pytest.param('N[C@@H](CCN1C(=O)C=CC1=O)C(=O)O', None, None,
                {4: 'maleimide_c'}, id='short_maleimide'),
            pytest.param('N[C@@H](CCCN1C(=O)C=CC1=O)C(=O)O', None, None,
                {4: 'maleimide_c'}, id='medium_maleimide'),
            pytest.param('N[C@H](CCCCN1C(=O)C=CC1=O)C(=O)O', None, None,
                {4: 'maleimide_c'}, id='d_lys_maleimide'),
            pytest.param('N[C@@H](CCCC)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='nle_no_sidechain'),
            pytest.param('N[C@@H](CCC)C(=O)O', {1: 7, 2: 6, 3: 7}, None, None,
                id='nva_no_sidechain'),
            pytest.param('CC(=O)Cl', {2: 6}, {2: '[Cl]'}, {2: 'backbone_c'},
                id='acetyl_chloride_cap'),
            pytest.param('O=C(Cl)c1ccccc1', {2: 6}, {2: '[Cl]'}, None,
                id='benzoyl_chloride_cap'),
            pytest.param('Cc1ccc(S(=O)(=O)Cl)cc1', {2: 16}, {2: '[Cl]'}, None,
                id='tosyl_chloride_cap'),
            pytest.param('CI', {2: 6}, {2: '[I]'}, None, id='methyl_iodide_cap'),
            pytest.param('CCBr', {2: 6}, {2: '[Br]'}, None, id='ethyl_bromide_cap'),
            pytest.param('ClCc1ccccc1', {2: 6}, {2: '[Cl]'}, None,
                id='benzyl_chloride_cap'),
            pytest.param('N', {1: 7}, {1: '[H]'}, {1: 'backbone_n'},
                id='ammonia_nucleophilic_cap'),
            pytest.param('CO', {1: 8}, {1: '[H]'}, None,
                id='methanol_nucleophilic_cap'),
            pytest.param('CC(=O)O', {2: 6}, {2: '[OH]'}, None,
                id='free_acid_leaving'),
        ],
    )
    def test_activation_slot_contract(
        self, smiles, expect_slots, expect_lg, expect_ct
    ):
        self._check(smiles, expect_slots, expect_lg, expect_ct)

    def test_asp_backbone_cooh_oh_excluded(self):
        """Asp: backbone COOH OH is not spuriously labelled as carboxyl."""
        r = self._check('N[C@@H](CC(=O)O)C(=O)O',
            expect_slots={1: 7, 2: 6, 3: 7, 4: 6}, expect_lg={4: '[OH]'},
            expect_ct={4: 'carboxyl', 1: 'backbone_n', 2: 'backbone_c'})
        carboxyl_slots = [s for s, ct in r.chem_types.items() if ct == 'carboxyl']
        assert len(carboxyl_slots) == 1, f"Expected 1 carboxyl slot, got {carboxyl_slots}"

    def test_glu_backbone_cooh_oh_excluded(self):
        """Glu: backbone COOH OH excluded from sidechain detection."""
        r = self._check('N[C@@H](CCC(=O)O)C(=O)O',
            expect_slots={1: 7, 2: 6, 3: 7, 4: 6}, expect_lg={4: '[OH]'},
            expect_ct={4: 'carboxyl', 1: 'backbone_n', 2: 'backbone_c'})
        carboxyl_slots = [s for s, ct in r.chem_types.items() if ct == 'carboxyl']
        assert len(carboxyl_slots) == 1

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

    def test_arg_guanidinium_slots(self):
        """All guanidine nitrogens retain their slots and actual functionality."""
        r = self._check('N[C@@H](CCCNC(=N)N)C(=O)O')
        assert r.chem_types == {1: 'backbone_n', 2: 'backbone_c',
                               3: 'backbone_n_mod', 4: 'guanidinium',
                               5: 'guanidinium', 6: 'guanidinium_imine',
                               7: 'guanidinium'}
        smap = self._slot_map(r.chuckles)
        guan_slot = [s for s, ct in r.chem_types.items() if ct == 'guanidinium'][0]
        assert smap[guan_slot] == 7
        assert r.leaving[guan_slot] == '[H]'

    def test_d_arg_guanidinium(self):
        """D-Arg has guanidine rather than sidechain primary amine sites."""
        r = self._check('N[C@H](CCCNC(=N)N)C(=O)O')
        assert 'guanidinium' in r.chem_types.values()
        assert 'amine_primary' not in r.chem_types.values()

    @pytest.mark.parametrize(
        'smiles,expected,r1_type',
        [
            ('O=C(O)[C@@H]1CCC(=O)N1',
             '[1*]N1C(=O)CC[C@H]1C([2*])=O', 'amide_nh'),
            ('O[C@@H](CC(N)=O)C(=O)O',
             '[1*]O[C@@H](CC(=O)N([4*])[5*])C([2*])=O', 'backbone_o'),
            ('N[C@@H](CCCNC(=N)N)C(=O)O',
             '[1*]N([3*])[C@@H](CCCN([5*])C(=N[6*])N([4*])[7*])C([2*])=O',
             'backbone_n'),
        ],
    )
    def test_backbone_roles_and_conjugated_n_preserve_selected_atoms(
            self, smiles, expected, r1_type):
        from pyPept.leaving_groups import restore_leaving_groups

        result = self._check(smiles)
        assert result.chuckles == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
        assert result.chem_types[1] == r1_type
        restored = restore_leaving_groups(Chem.MolFromSmiles(result.chuckles),
                                         result.leaving)
        assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

    def test_homoarg_guanidinium(self):
        """Homoarginine: longer chain to guanidinium."""
        r = self._check('N[C@@H](CCCCNC(=N)N)C(=O)O')
        assert 'guanidinium' in r.chem_types.values()

    def test_dopa_two_phenolic_oh(self):
        """DOPA (3,4-dihydroxyphenylalanine): two phenolic OHs detected."""
        r = self._check('N[C@@H](Cc1cc(O)c(O)cc1)C(=O)O')
        phenol_slots = [s for s, ct in r.chem_types.items() if ct == 'hydroxyl_phenolic']
        assert len(phenol_slots) == 2

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

    def test_phosphoserine_phosphate(self):
        """pSer: phosphate on serine OH."""
        r = self._check('N[C@@H](COP(=O)(O)O)C(=O)O')
        assert 'phosphate_p' in r.chem_types.values()
        smap = self._slot_map(r.chuckles)
        p_slot = [s for s, ct in r.chem_types.items() if ct == 'phosphate_p'][0]
        assert smap[p_slot] == 15
        assert r.leaving[p_slot] == '[OH]'

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

    def test_alkyne_dummy_on_carbon(self):
        """Alkyne: dummy lands on sp3 C (element 6)."""
        r = self._check('N[C@@H](CC#C)C(=O)O', expect_ct={4: 'alkyne_c'})
        smap = self._slot_map(r.chuckles)
        alk_slot = [s for s, ct in r.chem_types.items() if ct == 'alkyne_c'][0]
        assert smap[alk_slot] == 6

    def test_terminal_alkene_dummy_on_carbon(self):
        """Terminal alkene: dummy on internal vinyl C (element 6)."""
        r = self._check('N[C@@H](CCCC=C)C(=O)O', expect_ct={4: 'terminal_alkene'})
        smap = self._slot_map(r.chuckles)
        alk_slot = [s for s, ct in r.chem_types.items() if ct == 'terminal_alkene'][0]
        assert smap[alk_slot] == 6

    def test_tetrazine_dummy_on_carbon(self):
        """Tetrazine: dummy on sp3 C (element 6)."""
        r = self._check('N[C@@H](Cc1nncnn1)C(=O)O', expect_ct={4: 'tetrazine_c'})
        smap = self._slot_map(r.chuckles)
        tz_slot = [s for s, ct in r.chem_types.items() if ct == 'tetrazine_c'][0]
        assert smap[tz_slot] == 6

    def test_cyclooctyne_dummy_on_carbon(self):
        """Cyclooctyne: dummy on ring sp3 C (element 6)."""
        r = self._check('N[C@@H](CC1CCCCC#CC1)C(=O)O', expect_ct={4: 'cyclooctyne_c'})
        smap = self._slot_map(r.chuckles)
        cyo_slot = [s for s, ct in r.chem_types.items() if ct == 'cyclooctyne_c'][0]
        assert smap[cyo_slot] == 6

    def test_aldehyde_dummy_on_carbon(self):
        """Aldehyde: dummy on aldehyde C (element 6)."""
        r = self._check('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O', expect_ct={4: 'aldehyde'})
        smap = self._slot_map(r.chuckles)
        ald_slot = [s for s, ct in r.chem_types.items() if ct == 'aldehyde'][0]
        assert smap[ald_slot] == 6

    def test_nhs_ester_dummy_on_carbon(self):
        """NHS ester: dummy on sp3 C (element 6)."""
        r = self._check('N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O',
            expect_ct={4: 'nhs_ester'})
        smap = self._slot_map(r.chuckles)
        nhs_slot = [s for s, ct in r.chem_types.items() if ct == 'nhs_ester'][0]
        assert smap[nhs_slot] == 6

    def test_maleimide_dummy_on_carbon(self):
        """Maleimide: dummy on ring vinyl C (element 6)."""
        r = self._check('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O',
            expect_ct={4: 'maleimide_c'})
        smap = self._slot_map(r.chuckles)
        mal_slot = [s for s, ct in r.chem_types.items() if ct == 'maleimide_c'][0]
        assert smap[mal_slot] == 6

    @pytest.mark.parametrize(
        'smiles',
        [
            pytest.param('N[C@@H](CC1=CCCCCCC1)C(=O)O', id='tco_1c_linker'),
            pytest.param('N[C@@H](CCC1=CCCCCCC1)C(=O)O', id='tco_2c_linker'),
            pytest.param('N[C@@H](CCCCC1=CCCCCCC1)C(=O)O', id='tco_4c_linker'),
        ],
    )
    def test_tco_all_linkers_attach_outside_ring(self, smiles):
        result = self._check(smiles, expect_ct={4: 'tco_c'})
        molecule = Chem.MolFromSmiles(result.chuckles)
        dummy = next(atom for atom in molecule.GetAtoms()
                     if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 4)
        neighbor = dummy.GetNeighbors()[0]
        assert not neighbor.IsInRing(), "TCO dummy should be on exocyclic C"
