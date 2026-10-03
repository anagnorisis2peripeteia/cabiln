"""Attachment-site inference and reaction-family product chemistry."""

import os
import sys

import pytest
from hypothesis import given, strategies as st
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


from _chemistry_oracles import _romol
from _chemistry_fuzz import reordered
from _fuzzing import fuzz_settings, record


def test_explicit_trityl_protected_nitrogen_cannot_bypass_eligibility():
    from pyPept.attachments import attachment_site, reaction_for_types

    molecule = Chem.MolFromSmiles("CCN([4*])C(c1ccccc1)(c1ccccc1)c1ccccc1")
    site = attachment_site(molecule, 4)
    assert site["chem_type"] == "protected_amine"
    assert site["supported"] is False
    assert reaction_for_types(site["chem_type"], "nhs_ester") is None


def test_reaction_registry_rejects_collisions_and_resolves_aliases():
    import yaml
    from pyPept.interfaces.reaction_library import _YAML_PATH, compile_reactions

    entries = yaml.safe_load(_YAML_PATH.read_text())
    reactions, forward = compile_reactions(entries)
    _, reverse = compile_reactions(reversed(entries))
    assert {pair: entry['id'] for pair, entry in forward.items()} == {
        pair: entry['id'] for pair, entry in reverse.items()}
    assert reactions['aspartimide']['steps'] == reactions['imide_n_acylation']['steps']
    assert reactions['macrolactam_amide']['steps'] == reactions['backbone_amide']['steps']
    with pytest.raises(ValueError, match='Duplicate reaction route'):
        compile_reactions([*entries, {**entries[0], 'id': 'conflicting_amide'}])
    with pytest.raises(ValueError, match='duplicate reaction id'):
        compile_reactions([*entries, entries[0]])
    with pytest.raises(ValueError, match='cyclic reaction alias'):
        compile_reactions([{'id': 'a', 'alias_of': 'b'}, {'id': 'b', 'alias_of': 'a'}])


# Independent molecular expectations. {tail} is a spectator substituent, never
# generated from reaction SMARTS. Activated handles are explicitly supplied.
_FUZZ_PRODUCTS = [
    ("backbone_amide", "{tail}N[400*]", "[401*]C(=O)CC", "{tail}NC(=O)CC"),
    ("sidechain_amide", "{tail}N[400*]", "[401*]C(=O)CC", "{tail}NC(=O)CC"),
    ("isopeptide", "{tail}C(=O)[400*]", "[401*]NCC", "{tail}C(=O)NCC"),
    ("aspartimide", "{tail}C(=O)[400*]", "[401*]N(C)C(=O)C", "{tail}C(=O)N(C)C(=O)C"),
    ("hydroxylamine_cap", "{tail}O[400*]", "[401*]NCC", "{tail}ONCC"),
    ("sulfonamide", "{tail}N[400*]", "[401*]S(=O)(=O)CC", "{tail}NS(=O)(=O)CC"),
    ("backbone_ester", "{tail}O[400*]", "[401*]C(=O)CC", "{tail}OC(=O)CC"),
    ("disulfide", "{tail}S[400*]", "[401*]SCC", "{tail}SSCC"),
    ("diselenide", "{tail}[Se][400*]", "[401*][Se]CC", "{tail}[Se][Se]CC"),
    ("thioester", "{tail}S[400*]", "[401*]C(=O)CC", "{tail}SC(=O)CC"),
    ("thioester_sc", "{tail}S[400*]", "[401*]C(=O)CC", "{tail}SC(=O)CC"),
    ("thioether_halide", "{tail}S[400*]", "[401*]CCC", "{tail}SCCC"),
    ("thia_michael", "{tail}S[400*]", "[401*]CCC(=O)N", "{tail}SCCC(=O)N"),
    ("thiol_maleimide", "{tail}S[400*]", "[401*]C1=CC(=O)NC1=O", "{tail}SC1CC(=O)NC1=O"),
    ("nhs_ester_amide", "{tail}C([400*])C(=O)ON1C(=O)CCC1=O", "[401*]NCC", "{tail}CC(=O)NCC"),
    ("oxime_ligation", "{tail}O[NH][400*]", "CC([401*])=O", "{tail}ON=CC"),
    ("hydrazone", "{tail}N([400*])N", "CC([401*])=O", "{tail}NN=CC"),
    ("cuaac_1_4_triazole", "{tail}C([400*])C#C", "CC([401*])N=[N+]=[N-]", "{tail}Cc1cn(CC)nn1"),
    ("spaac_triazole", "{tail}C1([400*])C#CCCCCC1", "CC([401*])N=[N+]=[N-]", "{tail}C1c2nnn(CC)c2CCCCC1"),
    ("iedda_tetrazine_tco", "{tail}C([400*])c1nncnn1", "C1CC([401*])C=CCCC1", "{tail}CC1N=NC=C2CCCCCCC12"),
    ("phosphorylation", "{tail}O[400*]", "[401*]P(=O)(O)O", "{tail}OP(=O)(O)O"),
    ("rcm_alkene", "{tail}C([400*])=C", "CCC([401*])=C", "{tail}C=CCC"),
    ("imide_n_acylation", "{tail}N([400*])C(=O)C", "[401*]C(=O)CC", "{tail}N(C(=O)C)C(=O)CC"),
    ("guanidine_n_acylation", "{tail}NC(N[400*])=N", "[401*]C(=O)CC", "{tail}NC(NC(=O)CC)=N"),
    ("backbone_n_alkylation", "{tail}N([400*])C(=O)C", "[401*]CCC", "{tail}N(CCC)C(=O)C"),
    ("macrolactam_amide", "{tail}N[400*]", "[401*]C(=O)CC", "{tail}NC(=O)CC"),
    ("ester_sidechain_oh", "{tail}O[400*]", "[401*]C(=O)CC", "{tail}OC(=O)CC"),
    ("n_alkylation_halide_backbone", "{tail}C[400*]", "[401*]NCC", "{tail}CNCC"),
    ("n_o_bond", "{tail}N[400*]", "[401*]OCC", "{tail}NOCC"),
    ("n_n_bond", "{tail}N[400*]", "[401*]NCC", "{tail}NNCC"),
    ("aryl_c_c_bond", "{tail}C[400*]", "[401*]c1ccccc1", "{tail}Cc1ccccc1"),
    ("aryl_amide_to_backbone_n", "{tail}c1ccc(C(=O)[400*])cc1", "[401*]NCC", "{tail}c1ccc(C(=O)NCC)cc1"),
    ("reduced_amide", "{tail}N[400*]", "[401*]CCC", "{tail}NCCC"),
    ("benzylamine_xlink", "{tail}c1ccc(C[400*])cc1", "[401*]NCC", "{tail}c1ccc(CNCC)cc1"),
    ("aryl_ether", "{tail}c1ccc(O[400*])cc1", "[401*]C(C)(C)CC", "{tail}c1ccc(OC(C)(C)CC)cc1"),
    ("aryl_o_alkylation", "{tail}c1ccc(O[400*])cc1", "[401*]CCC", "{tail}c1ccc(OCCC)cc1"),
    ("aryl_ester", "{tail}c1ccc(O[400*])cc1", "[401*]C(=O)CC", "{tail}c1ccc(OC(=O)CC)cc1"),
    ("diaryl_ether", "{tail}c1ccc(O[400*])cc1", "[401*]c1ccccc1", "{tail}c1ccc(Oc2ccccc2)cc1"),
]


@pytest.mark.fuzz
@pytest.mark.parametrize("reaction_id,left,right,expected", _FUZZ_PRODUCTS,
                         ids=[row[0] for row in _FUZZ_PRODUCTS])
@fuzz_settings(examples=10)
@given(tail=st.sampled_from(("C", "CCC", "[13CH3]", "[NH3+]C", "F[C@H](C)", "F[C@@H](C)")),
       order=st.integers(0, 65535), port=st.integers(0, 100), reverse=st.booleans())
def test_fuzz_reaction_products_preserve_spectators_and_target_only_selected_ports(
    reaction_id, left, right, expected, tail, order, port, reverse
):
    from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks

    record("reaction." + reaction_id, tail=tail, order=order, port=port, reverse=reverse)
    molecules = []
    first, second = 400 + port, 800 + port
    for source, old, new in ((left, 400, first), (right, 401, second)):
        text = source.format(tail=tail).replace(f"[{old}*]", f"[{new}*]")
        molecule = Chem.MolFromSmiles(text)
        assert molecule is not None, text
        molecules.append(reordered(molecule, order))
    if reverse:
        molecules.reverse()
        first, second = second, first
    product = run_bond_smirks(
        *molecules, first, second, REACTIONS[reaction_id], False, connection_id=7
    )
    literal = Chem.MolFromSmiles(expected.format(tail=tail))
    assert literal is not None
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(literal)
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())
    assert any(
        bond.HasProp("_connection_idx") and bond.GetIntProp("_connection_idx") == 7
        for bond in product.GetBonds()
    )
    # A requested port absent from the molecule must not silently use another.
    with pytest.raises(ValueError, match="produced no products"):
        run_bond_smirks(*molecules, first + 2000, second, REACTIONS[reaction_id], False)


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


class TestActivatedSiteInference:
    """Infer site chemistry from already activated fragments."""

    @pytest.mark.parametrize(
        'smiles,expected',
        [
            pytest.param('C([4*])=O', 'aldehyde', id='aldehyde-one-hydrogen'),
            pytest.param('CC([4*])=O', 'aldehyde', id='aldehyde-no-hydrogen'),
            pytest.param('[4*][NH]OC', 'aminooxy', id='aminooxy'),
            pytest.param('[4*]N(N)C=O', 'hydrazide', id='hydrazide'),
            pytest.param('C([4*])1=CC(=O)NC1=O', 'maleimide_c', id='maleimide'),
            pytest.param('C([4*])C(=O)ON1C(=O)CCC1=O', 'nhs_ester', id='nhs-ester'),
            pytest.param('c1ccc2n([4*])ccc2c1', 'aromatic_nh', id='indole-n'),
            pytest.param('CC(=O)N([4*])[5*]', 'amide_nh', id='primary-amide-two-sites'),
            pytest.param('CC(=O)N([4*])C', 'amide_nh', id='secondary-amide'),
            pytest.param('CS(=O)(=O)N([4*])C', 'sulfonamide_nh', id='sulfonamide'),
            pytest.param('COC(=O)N([4*])C', 'carbamate_nh', id='carbamate'),
            pytest.param('CNC(=O)N([4*])C', 'urea_nh', id='urea'),
            pytest.param('COP(=O)(O)N([4*])C', 'phosphoramide_nh', id='phosphoramide'),
            pytest.param('CC(=N)N([4*])C', 'amidine_nh', id='amidine-amino-n'),
            pytest.param('[4*]N=C(N)C', 'amidine_imine', id='amidine-imine-n'),
            pytest.param('CS(=O)(=O)[4*]', 'sulfonyl_s', id='sulfonyl-donor'),
            pytest.param('CON([4*])C', 'substituted_aminooxy', id='n-methyl-aminooxy'),
            pytest.param('CO[N+]([4*])(C)C', 'element_7', id='charged-aminooxy'),
            pytest.param('C[NH2+][4*]', 'element_7', id='ammonium-not-neutral-amine'),
            pytest.param('[4*][NH+]([5*])C', 'element_7', id='two-slot-ammonium'),
            pytest.param('O=P([4*])([5*])[6*]', 'phosphate_p', id='three-phosphate-ports'),
            pytest.param('NC(=N)N([4*])C', 'guanidinium', id='guanidine-chain-n'),
            pytest.param('N=C(N([4*])[5*])NC', 'guanidinium', id='guanidine-terminal-n'),
            pytest.param('COP([4*])(=O)O', 'phosphate_p', id='activated-phosphate'),
            pytest.param('CN([5*])C([4*])=O', 'formamide_c', id='formamide-not-aldehyde'),
            pytest.param('C[N+]([4*])(C)C', 'element_7', id='quaternary-n-not-primary-amine'),
            pytest.param('[4*]N=CC', 'element_7', id='imine-not-primary-amine'),
        ],
    )
    def test_activated_site_type(self, smiles, expected):
        from pyPept.interfaces.reaction_library import infer_chem_type

        molecule = Chem.MolFromSmiles(smiles)
        assert molecule is not None
        dummy = next(atom for atom in molecule.GetAtoms()
                     if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 4)
        attachment = dummy.GetNeighbors()[0].GetIdx()
        assert infer_chem_type(molecule, attachment, slot=3) == expected

    @pytest.mark.parametrize('smiles,expected', [
        ('CN([4*])C(=O)N(C)[5*]', 'urea_nh'),
        ('CN([4*])C(=N)N(C)[5*]', 'guanidinium'),
    ])
    def test_symmetric_functional_groups_classify_every_anchor(self, smiles, expected):
        from pyPept.site_chemistry import Perception

        molecule = Chem.MolFromSmiles(smiles)
        full = Perception(molecule)
        targeted = Perception(Chem.Mol(molecule), targeted=True)
        for dummy in molecule.GetAtoms():
            if dummy.GetAtomicNum() == 0:
                index = dummy.GetNeighbors()[0].GetIdx()
                assert full.nitrogen(index) == targeted.nitrogen(index) == expected

    @pytest.mark.parametrize(
        'template,slot,declaration,leaving,expected',
        [
            ('[1*]Cc1ccc(C[2*])cc1', 1, '1:backbone_c_red', '[H]', 'backbone_c_red'),
            ('[1*]N([3*])[C@@H]([2*])C(C)C', 2, '2:sp3_c_anchor', '[H]', 'sp3_c_anchor'),
            ('[4*]C(=O)c1ccccc1', 4, '', '[OH]', 'aryl_amide_c'),
            ('[4*]N1C=CC=C1', 4, '4:amine_primary', '[H]', 'aromatic_nh'),
            ('CC(=O)N([4*])[5*]', 4, '4:backbone_n', '[H]', 'amide_nh'),
            ('[4*]OC', 4, '4:backbone_c_red', '[H]', 'hydroxyl'),
        ],
    )
    def test_structurally_checked_declarations(self, template, slot, declaration,
                                               leaving, expected):
        from pyPept.attachments import attachment_site

        molecule = Chem.MolFromSmiles(template)
        molecule.SetProp('m_chem_types', declaration)
        groups = [None] * slot
        groups[slot - 1] = leaving
        assert attachment_site(molecule, slot, groups)['chem_type'] == expected

    @pytest.mark.parametrize(
        'notation,expected',
        [
            ('N.Ac(4,2)', 'CC(=O)NC(=O)C[C@H](N)C(=O)O'),
            ('R.Ac(7,2)', 'CC(=O)N(CCC[C@H](N)C(=O)O)C(=N)N'),
            ('N.D(4,4)', 'N[C@@H](CC(=O)NC(=O)C[C@H](N)C(=O)O)C(=O)O'),
            ('R.D(7,4)', 'N=C(N)N(CCC[C@H](N)C(=O)O)C(=O)C[C@H](N)C(=O)O'),
            ('ac-Pyr-am', 'CC(=O)N1C(=O)CC[C@H]1C(N)=O'),
            ('K.Ac(5,2)', 'CC(=O)NCCCC[C@H](N)C(=O)O'),
            ('pXyl.A(1,1)', 'Cc1ccc(CN[C@@H](C)C(=O)O)cc1'),
            ('pXyl.A(2,1)', 'Cc1ccc(CN[C@@H](C)C(=O)O)cc1'),
            ('ValAryl.ImzScaffold(2,5)', 'N[C@@H](c1ncn(C)c1)C(C)C'),
            ('Ser_PO3H2.S(4,4)', 'N[C@@H](COP(=O)(O)OC[C@H](N)C(=O)O)C(=O)O'),
        ],
    )
    def test_corrected_types_preserve_explicit_products(self, notation, expected):
        assert Chem.MolToSmiles(_romol(notation)) == Chem.MolToSmiles(
            Chem.MolFromSmiles(expected))

    @pytest.mark.parametrize('notation', ['W.Ac(4,2)', 'H.Ac(4,2)'])
    def test_aromatic_n_does_not_advertise_aliphatic_amine_reaction(self, notation):
        with pytest.raises(ValueError, match='aromatic_nh') as error:
            _romol(notation)
        assert 'produced no products' not in str(error.value)

    @pytest.mark.parametrize(
        'template,slot,declaration,compatible',
        [
            ('[4*]Oc1ccccc1', 4, '4:hydroxyl_phenolic', True),
            ('CN([4*])[5*]', 5, '5:amine_secondary', True),
            ('CN([4*])[5*]', 4, '4:amine_secondary', False),
            ('CC(=O)N([4*])[5*]', 5, '5:amine_secondary', False),
            ('N=C(N([4*])[5*])NC', 5, '5:amine_secondary', False),
            ('CCO[2*]', 2, '2:backbone_o', True),
            ('CCO[4*]', 4, '4:backbone_o', False),
        ],
    )
    def test_legacy_role_aliases_require_their_structural_context(
            self, template, slot, declaration, compatible):
        from pyPept.attachments import attachment_site, declaration_is_compatible

        molecule = Chem.MolFromSmiles(template)
        molecule.SetProp('m_chem_types', declaration)
        site = attachment_site(molecule, slot, ['[H]'] * slot)
        assert declaration_is_compatible(molecule, site) is compatible

    @pytest.mark.parametrize(
        'template,leaving,expected',
        [('C[4*]', '[Cl]', 'CSC'),
         ('C([4*])Cl', '[Cl]', 'CSCCl'),
         ('C([4*])Cl', '[H]', None),
         ('C([4*])Cl', None, None)],
    )
    def test_halide_substitution_uses_selected_leaving_group(
            self, template, leaving, expected):
        from pyPept.attachments import resolve_connection
        from pyPept.interfaces.reaction_library import run_bond_smirks

        thiol = Chem.MolFromSmiles('CS[3*]')
        electrophile = Chem.MolFromSmiles(template)
        electrophile.SetProp('m_chem_types', '4:alkyl_halide_c')
        connection = resolve_connection(
            thiol, 3, electrophile, 4,
            leaving_groups1=[None, None, '[H]'],
            leaving_groups2=[None, None, None, leaving])
        if expected is None:
            assert connection.chem_type2 != 'alkyl_halide_c'
            assert connection.reaction is None
        else:
            assert connection.chem_type2 == 'alkyl_halide_c'
            product = run_bond_smirks(thiol, electrophile, 3, 4,
                                     connection.reaction, intramolecular=False)
            assert Chem.MolToSmiles(product) == Chem.MolToSmiles(
                Chem.MolFromSmiles(expected))


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
        """One chloride is selected for substitution on a dichloromethyl group."""
        from pyPept.interfaces.reaction_library import infer_chem_type
        # CH2Cl2 -> [400*]CH2Cl, with the displaced chloride in metadata.
        mol = Chem.MolFromSmiles('[400*]CCl')
        assert mol is not None
        attach_idx = next(
            nb.GetIdx()
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 400
            for nb in a.GetNeighbors()
        )
        ct = infer_chem_type(mol, attach_idx, slot=400, leaving='[Cl]')
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

    def test_tz_tco_iedda_reaction_index(self):
        """REACTION_INDEX routes (tetrazine_c, tco_c) → iedda_tetrazine_tco (routing check)."""
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get(('tetrazine_c', 'tco_c'))
        assert entry is not None
        assert entry['id'] == 'iedda_tetrazine_tco'
        assert entry.get('take_largest') is True, "IEDDA must have take_largest=True (N2 byproduct)"


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
        'amine_secondary':  'N[C@@H](CCCCNC)C(=O)O',
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
