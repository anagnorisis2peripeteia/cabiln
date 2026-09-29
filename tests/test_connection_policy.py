"""Numbered-site eligibility agrees across parsing, the builder and assembly."""

import warnings

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import attachments, monomer_store
from pyPept.molecule import Molecule
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
    real_infer = attachments.infer_chem_type
    calls = []

    def infer(mol, index, **kwargs):
        calls.append((mol, kwargs["slot"]))
        return real_infer(mol, index, **kwargs)

    monkeypatch.setattr(attachments, "infer_chem_type", infer)
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


def test_supported_exotic_warning_still_reaches_validation_report():
    report = Sequence.validate("K.!n(4,4)%K.!n(4,4)")
    assert report.ok, report.errors
    assert len(report.warnings) == 1
    assert "N–N join (hydrazide/hydrazone)" in report.warnings[0]
    assert "R4 of residue 0 <-> R4 of residue 1" in report.warnings[0]
    assert Molecule(report.build()).get_molecule(fmt="ROMol") is not None


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
    with pytest.raises(ValueError, match="No reaction defined.*wrong R-group index"):
        Molecule(sequence)
