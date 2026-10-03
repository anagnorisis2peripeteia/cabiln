"""Numbered-site eligibility agrees across parsing, the builder and assembly."""

import warnings

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import attachments, monomer_store
from pyPept.molecule import Molecule
from pyPept.peptide import Endpoint
from pyPept.sequence import Sequence, _check_bond_chemistry
from pyPept.web.app import create_app


@pytest.fixture
def connection_library(tmp_path, monkeypatch):
    path = tmp_path / "connection-sites.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for symbol, smiles, slot, leaving in (
            ("HydroxyCarrier", "CCO[4*]", 4, "[H]"),
            ("PhosphateUnit", "[4*]P(=O)(O)O", 4, "[OH]"),
            ("HighAmine", "CN[165*]", 165, "[H]"),
            ("IsotopicAmine", "[1CH3]N[165*]", 165, "[H]"),
            ("LowAcid", "CC(=O)[65*]", 65, "[OH]"),
        ):
            mol = Chem.MolFromSmiles(smiles)
            for key, value in {
                "symbol": symbol,
                "m_abbr": symbol,
                "m_type": "chem",
                "m_subtype": "undefined",
                "m_Rgroups": ",".join(["None"] * (slot - 1) + [leaving]),
                # Descriptive metadata must not determine compatibility.
                "m_chem_types": f"{slot}:declared_label",
            }.items():
                mol.SetProp(key, value)
            writer.write(mol)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    yield
    monomer_store._invalidate_sdf()


@pytest.mark.parametrize("reverse", [False, True])
def test_phosphorylation_from_http_sites_through_notation_to_product(
    connection_library, reverse
):
    symbols = ["HydroxyCarrier", "PhosphateUnit"]
    if reverse:
        symbols.reverse()
    with TestClient(create_app()) as client:
        sites = []
        for symbol in symbols:
            response = client.get("/monomer_rgroups", params={"abbr": symbol})
            assert response.status_code == 200, response.text
            (site,) = response.json()["rgroups"]
            assert site["declared_chem_type"] == "declared_label"
            sites.append(site)
        assert {site["chem_type"] for site in sites} == {"hydroxyl", "phosphate_p"}
        response = client.post(
            "/validate_bond",
            json={
                "chem_type_a": sites[0]["chem_type"],
                "chem_type_b": sites[1]["chem_type"],
            },
        )
        assert response.status_code == 200, response.text
        assert response.json()["valid"] is True
        assert response.json()["reaction"] == "phosphorylation"
        inserted = client.post(
            "/insert_bond",
            json={
                "cabiln": symbols[0],
                "host_residue_idx": 0,
                "new_abbr": symbols[1],
                "r_host": 4,
                "r_new": 4,
            },
        )
        assert inserted.status_code == 200, inserted.text
    explicit = "%".join(f"{symbol}.!p(4,4)" for symbol in symbols)
    expected = Chem.MolToSmiles(Chem.MolFromSmiles("CCOP(=O)(O)O"))
    for source in (explicit, inserted.json()["result"]):
        with warnings.catch_warnings(record=True) as caught:
            sequence = Sequence(source)
        assert not caught
        product = Molecule(sequence).get_molecule(fmt="ROMol")
        assert Chem.MolToSmiles(product) == expected
        assert all(atom.GetAtomicNum() for atom in product.GetAtoms())


def test_resolve_connection_uses_only_selected_sites_and_current_leaving(monkeypatch):
    carbon = Chem.MolFromSmiles("[1*]NCC(=O)[4*]")
    carbon.SetProp("m_Rgroups", "[H],None,None,[H]")
    carbon.SetProp("m_chem_types", "1:declared_n,4:declared_c")
    nitrogen = Chem.MolFromSmiles("[1*]NC")
    real_classify = attachments.Perception.classify
    calls = []

    def classify(perception, index, slot=None, leaving=None, declared=None):
        calls.append((perception.mol, slot))
        return real_classify(perception, index, slot, leaving, declared)

    monkeypatch.setattr(attachments.Perception, "classify", classify)
    connection = attachments.resolve_connection(
        carbon,
        4,
        nitrogen,
        1,
        leaving_groups1=["[H]", None, None, "[OH]"],
        leaving_groups2="[H]",
    )
    assert calls == [(carbon, 4), (nitrogen, 1)]
    assert connection.chem_type1 == "carboxyl"
    assert connection.chem_type2 == "backbone_n"
    assert connection.reaction["id"] == "isopeptide"
    # The stored H still resolves as aldehyde; no caller override was cached.
    assert attachments.attachment_site(carbon, 4)["chem_type"] == "aldehyde"
    assert attachments.attachment_site(carbon, 4)["declared_chem_type"] == "declared_c"


def test_slots_on_the_same_atom_keep_distinct_chemistry():
    nitrogen = Chem.MolFromSmiles("[1*]N([3*])C")
    carbon = Chem.MolFromSmiles("CC(=O)[4*]")
    carbon.SetProp("m_Rgroups", "None,None,None,[OH]")
    backbone = attachments.resolve_connection(nitrogen, 1, carbon, 4)
    modification = attachments.resolve_connection(nitrogen, 3, carbon, 4)
    assert backbone.chem_type1 == "backbone_n"
    assert backbone.reaction["id"] == "isopeptide"
    assert modification.chem_type1 == "backbone_n_mod"
    assert modification.reaction["id"] == "imide_n_acylation"


def test_wrong_slot_diagnostics_and_http_rejection(connection_library):
    missing = Sequence.validate("HydroxyCarrier.!p(7,4)%PhosphateUnit.!p(4,7)")
    assert not missing.ok
    assert "Residue 0 has no R7 attachment point" in missing.errors[0]
    unsupported = Sequence.validate("HydroxyCarrier.!p(4,4)%HydroxyCarrier.!p(4,4)")
    assert not unsupported.ok
    assert "O–O inter-monomer bond" in unsupported.errors[0]
    assert "correct R-group numbers" in unsupported.errors[0]
    with TestClient(create_app()) as client:
        response = client.post(
            "/validate_bond",
            json={
                "chem_type_a": "hydroxyl",
                "chem_type_b": "hydroxyl",
            },
        )
    assert response.json() == {
        "valid": False,
        "reason": "No reaction for hydroxyl + hydroxyl",
    }


@pytest.mark.parametrize(
    "source, expected, ownership",
    [
        ("HighAmine%LowAcid", "CN.CC(=O)O", [2, 4]),
        ("LowAcid%HighAmine", "CN.CC(=O)O", [4, 2]),
        ("HighAmine.!a(165,65)%LowAcid.!a(65,165)", "CNC(C)=O", [2, 3]),
        # A temporary dummy label must not consume an isotope-labeled real atom.
        ("IsotopicAmine.!a(165,65)%LowAcid.!a(65,165)", "[1CH3]NC(C)=O", [2, 3]),
    ],
)
def test_large_slots_keep_distinct_endpoints_and_restoration(
    connection_library, source, expected, ownership
):
    assembled = Molecule(Sequence(source), depiction=None)
    product = assembled.get_molecule(fmt="ROMol")
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    owners = assembled.get_residue_atom_map()
    assert [len(owners[index]) for index in range(2)] == ownership
    assert sorted(atom for indices in owners.values() for atom in indices) == list(
        range(product.GetNumAtoms())
    )


def test_bare_molecule_compatibility_does_not_claim_numbered_site_support():
    oxygen = Chem.MolFromSmiles("CO")
    phosphorus = Chem.MolFromSmiles("OP(=O)(O)O")
    with pytest.raises(ValueError, match="O–15 inter-monomer bond"):
        _check_bond_chemistry(oxygen, 1, phosphorus, 1)


def test_supported_numbered_reaction_does_not_create_an_atom_pair_warning():
    report = Sequence.validate("K.!n(4,4)%K.!n(4,4)")
    assert report.ok, report.errors
    assert report.warnings == []
    product = Molecule(report.build()).get_molecule(fmt="ROMol")
    assert len(product.GetSubstructMatches(Chem.MolFromSmarts("[N]-[N]"))) == 1
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())


@pytest.mark.parametrize("tag", ["1", "n"])
def test_legacy_parseable_graph_does_not_claim_reaction_support(tag):
    report = Sequence.validate(f"K.!{tag}(3,3)-A-K.!{tag}(3,3)-am")
    assert report.ok, report.errors
    assert any("N–N join" in warning for warning in report.warnings)
    sequence = report.build()
    assert sequence.s_nmonomers == 4
    connection = attachments.resolve_connection(
        sequence.s_monomers[0]["m_romol"],
        3,
        sequence.s_monomers[2]["m_romol"],
        3,
    )
    assert connection.chem_type1 == connection.chem_type2 == "backbone_n_mod"
    assert connection.reaction is None
    with TestClient(create_app()) as client:
        response = client.post(
            "/validate_bond",
            json={
                "chem_type_a": connection.chem_type1,
                "chem_type_b": connection.chem_type2,
            },
        )
    assert response.json()["valid"] is False
    with pytest.raises(ValueError, match="No reaction defined.*backbone_n_mod.*slot 3.*not supported"):
        Molecule(sequence)


@pytest.fixture
def substitution_library(tmp_path, monkeypatch):
    path = tmp_path / "substitution.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for symbol, smiles, groups in (
            ("Amine", "CN([4*])[5*]", "None,None,None,[H],[H]"),
            ("Methyl", "C[4*]", "None,None,None,[Br]"),
            ("Acyl", "CC(=O)[4*]", "None,None,None,[OH]"),
            ("NHS", "C([4*])C(=O)ON1C(=O)CCC1=O", "None,None,None,[H]"),
            ("Alkyne", "C([4*])([5*])C#C", "None,None,None,[H],[H]"),
            ("Azide", "C([4*])N=[N+]=[N-]", "None,None,None,[H]"),
            ("Sulfonyl", "CS(=O)(=O)[4*]", "None,None,None,[OH]"),
            ("Urea", "CN([6*])C(=O)N([4*])[5*]", "None,None,None,[H],[H],[H]"),
            ("Carbamate", "COC(=O)N([4*])[5*]", "None,None,None,[H],[H]"),
            ("Hydrazide", "CN([5*])CC(=O)N([4*])N", "None,None,None,[H],[H]"),
            ("Aminooxy", "CN([5*])CCON[4*]", "None,None,None,[H],[H]"),
            ("Aldehyde", "CC(=O)[4*]", "None,None,None,[H]"),
            ("Phosphate", "O=P([4*])([5*])[6*]", "None,None,None,[OH],[OH],[OH]"),
            ("ChargedPhosphate", "O=P([4*])([5*])[O-]", "None,None,None,[OH],[OH]"),
            ("Alcohol", "CCO[4*]", "None,None,None,[H]"),
            ("Guanidine", "N=C(N([4*])[5*])N(C)[6*]", "None,None,None,[H],[H],[H]"),
        ):
            molecule = Chem.MolFromSmiles(smiles)
            for key, value in {"symbol": symbol, "m_abbr": symbol,
                               "m_type": "chem", "m_subtype": "undefined",
                               "m_Rgroups": groups, "m_chem_types": {
                                   "Amine": "4:amine_primary,5:amine_primary",
                                   "Methyl": "4:alkyl_halide_c", "Acyl": "4:carboxyl",
                                   "NHS": "4:nhs_ester", "Alkyne": "4:alkyne_c,5:alkyne_c",
                                   "Azide": "4:azide_alpha_c",
                               }.get(symbol, "")}.items():
                molecule.SetProp(key, value)
            writer.write(molecule)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    yield
    monomer_store._invalidate_sdf()


@pytest.mark.parametrize("source,kind,used,expected", [
    ("Amine", "amine_primary", [False, False], "CN"),
    ("Amine.Methyl(4,4)", "amine_secondary", [True, False], "CNC"),
    ("Amine.Acyl(4,4)", "amide_nh", [True, False], "CNC(C)=O"),
    ("Amine.Methyl(4,4).Methyl(5,4)", "element_7", [True, True], "CN(C)C"),
])
def test_current_chemistry_keeps_numbered_slots_after_substitution(
    substitution_library, source, kind, used, expected
):
    assembly = Molecule(Sequence(source), depiction=None)
    sites = assembly.current_sites(0)
    assert [site["slot"] for site in sites] == [4, 5]
    assert [site["used"] for site in sites] == used
    assert [site["functionality"] for site in sites] == [kind, kind]
    assert Chem.MolToSmiles(assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    anchors = assembly.get_attachment_atom_map()
    assert anchors[Endpoint(0, 4)] == anchors[Endpoint(0, 5)]


@pytest.mark.parametrize('source,kind', [
    ('Amine.Acyl(4,4)', 'amide_nh'),
    ('Amine.Sulfonyl(4,4)', 'sulfonamide_nh'),
])
def test_builder_checks_current_structure_instead_of_submitted_type_labels(
    substitution_library, source, kind
):
    with TestClient(create_app()) as client:
        sites = client.get("/monomer_rgroups", params={
            "abbr": "Amine", "residue_idx": 0, "cabiln": source,
        }).json()["rgroups"]
        assert sites[1]["chem_type"] == kind
        result = client.post("/validate_bond", json={
            "cabiln": source, "residue_idx_a": 0, "slot_a": 4,
            "abbr_b": "Methyl", "slot_b": 4,
            "chem_type_a": "amine_primary", "chem_type_b": "alkyl_halide_c",
        })
        assert result.status_code == 200
        assert result.json()["valid"] is False
        assert "already" in result.json()["reason"]
        request = {"cabiln": source, "residue_idx_a": 0, "slot_a": 5,
                   "abbr_b": "NHS", "slot_b": 4,
                   "chem_type_a": "amine_primary", "chem_type_b": "nhs_ester"}
        assert client.post('/validate_bond', json=request).json()['valid'] is False
        inserted = client.post('/insert_bond', json={
            'cabiln': source, 'host_residue_idx': 0, 'r_host': 5,
            'new_abbr': 'NHS', 'r_new': 4,
        })
        assert inserted.status_code == 400
        assert kind in inserted.json()['error']
        request['cabiln'] = 'Amine.Methyl(4,4)'
        assert client.post('/validate_bond', json=request).json()['valid'] is True


def test_builder_can_use_the_second_amine_slot_after_alkylation(substitution_library):
    from pyPept.editor import PeptideDocument

    first = PeptideDocument("Amine")
    source = first.attach(first.select(0), 4, "Methyl", 4)
    second = PeptideDocument(source)
    source = second.attach(second.select(0), 5, "Methyl", 4)
    assert Chem.MolToSmiles(Molecule(Sequence(source), depiction=None).mol) == "CN(C)C"


def test_shared_reactive_handle_is_consumed_even_when_another_slot_remains(substitution_library):
    from pyPept.editor import PeptideDocument

    before = PeptideDocument('Alkyne')
    assert before.current_site(before.select(0), 5)['supported']
    source = before.attach(before.select(0), 4, 'Azide', 4)
    after = PeptideDocument(source)
    remaining = after.current_site(after.select(0), 5)
    assert remaining['used'] is False
    assert remaining['supported'] is False
    assert 'consumed' in remaining['reason']
    assert len(after.assembly.mol.GetSubstructMatches(Chem.MolFromSmarts('n1nncc1'))) == 1
    with pytest.raises(ValueError, match='consumed'):
        after.attach(after.select(0), 5, 'Azide', 4)


def test_large_graph_checks_preserve_full_context_beyond_a_thousand_matches():
    from pyPept.site_chemistry import Perception

    molecule = Chem.MolFromSmiles('.'.join(['CNC(C)=O'] * 1100 + ['CNC', 'CN']))
    full = Perception(molecule)
    targeted = Perception(Chem.Mol(molecule), targeted=True)
    assert len(full.matches('amide_nh')) == 1100
    nitrogens = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 7]
    # The final amides used to fall past RDKit's default match limit.
    for position in (0, 999, 1000, 1099, 1100, 1101):
        expected = 'amide_nh' if position < 1100 else 'amine_secondary' if position == 1100 else 'amine_primary'
        assert full.nitrogen(nitrogens[position]) == expected
        assert targeted.nitrogen(nitrogens[position]) == expected


@pytest.mark.parametrize('source,kind,expected', [
    ('Amine.Sulfonyl(4,4)', 'sulfonamide_nh', 'CNS(C)(=O)=O'),
    ('Urea', 'urea_nh', 'CNC(N)=O'),
    ('Carbamate', 'carbamate_nh', 'COC(N)=O'),
])
def test_special_nitrogens_keep_alkylation_without_free_amine_coupling(
    substitution_library, source, kind, expected
):
    from pyPept.editor import PeptideDocument

    document = PeptideDocument(source)
    selection = document.select(0)
    site = document.current_site(selection, 5)
    assert site['functionality'] == site['chem_type'] == kind
    assert not site['used']
    assert Chem.MolToSmiles(document.assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    with pytest.raises(ValueError, match=kind):
        document.attach(selection, 5, 'NHS', 4)
    attached = PeptideDocument(document.attach(selection, 5, 'Methyl', 4))
    expected_product = {
        'sulfonamide_nh': 'CN(C)S(C)(=O)=O',
        'urea_nh': 'CNC(=O)NC',
        'carbamate_nh': 'COC(=O)NC',
    }[kind]
    assert Chem.MolToSmiles(attached.assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected_product))
    assert attached.current_site(attached.select(0), 5)['used']


@pytest.mark.parametrize('source,kind,expected', [
    ('C._Me(4,2)', 'thioether', 'CSC[C@H](N)C(=O)O'),
    ('C.ac(4,2)', 'thioester', 'CC(=O)SC[C@H](N)C(=O)O'),
    ('C.C(4,4)', 'disulfide', 'N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O'),
    ('Sec.Sec(4,4)', 'diselenide', 'N[C@@H](C[Se][Se]C[C@H](N)C(=O)O)C(=O)O'),
])
def test_used_sulfur_and_selenium_report_product_functionality(source, kind, expected):
    assembly = Molecule(Sequence(source), depiction=None)
    site = next(site for site in assembly.current_sites(0) if site['slot'] == 4)
    assert site['used']
    assert site['functionality'] == kind
    assert Chem.MolToSmiles(assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))


@pytest.mark.parametrize('symbol,expected', [
    ('Hydrazide', 'CNCC(=O)NN=CC'),
    ('Aminooxy', 'CNCCON=CC'),
])
def test_condensation_consumes_its_handle_and_preserves_a_separate_amine(
    substitution_library, symbol, expected
):
    from pyPept.editor import PeptideDocument

    before = PeptideDocument(symbol)
    after = PeptideDocument(before.attach(before.select(0), 4, 'Aldehyde', 4))
    assert Chem.MolToSmiles(after.assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    assert after.current_site(after.select(0), 4)['used']
    assert after.current_site(after.select(0), 5)['chem_type'] == 'amine_secondary'
    with pytest.raises(ValueError, match='already bonded'):
        after.attach(after.select(0), 4, 'Aldehyde', 4)
    assert after.check_connection(after.select(0), 5, 'NHS', 4)['id'] == 'nhs_ester_amide'


@pytest.mark.parametrize('symbol,slots,expected', [
    ('Phosphate', (4, 5, 6), 'CCOP(=O)(OCC)OCC'),
    ('ChargedPhosphate', (4, 5), 'CCOP(=O)([O-])OCC'),
])
def test_each_phosphate_port_is_independent_and_preserves_charge(
    substitution_library, symbol, slots, expected
):
    from pyPept.editor import PeptideDocument

    document = PeptideDocument(symbol)
    for index, slot in enumerate(slots):
        sites = document.assembly.current_sites(0)
        assert [site['slot'] for site in sites] == list(slots)
        assert [site['used'] for site in sites] == [position < index for position in range(len(slots))]
        assert all(site['chem_type'] == 'phosphate_p' and site['supported']
                   for site in sites if not site['used'])
        document = PeptideDocument(document.attach(document.select(0), slot, 'Alcohol', 4))
    assert Chem.MolToSmiles(document.assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    assert all(site['used'] for site in document.assembly.current_sites(0))


def test_guanidine_acylation_reclassifies_shared_nitrogen_only(substitution_library):
    from pyPept.editor import PeptideDocument

    before = PeptideDocument('Guanidine')
    after = PeptideDocument(before.attach(before.select(0), 4, 'Acyl', 4))
    assert after.current_site(after.select(0), 5)['chem_type'] == 'amide_nh'
    assert after.current_site(after.select(0), 6)['chem_type'] == 'guanidinium'
    with pytest.raises(ValueError, match='amide_nh'):
        after.check_connection(after.select(0), 5, 'NHS', 4)
    assert after.check_connection(after.select(0), 6, 'Acyl', 4)['id'] == 'guanidine_n_acylation'
    assert Chem.MolToSmiles(after.assembly.mol) == Chem.MolToSmiles(Chem.MolFromSmiles('CNC(=N)NC(C)=O'))
