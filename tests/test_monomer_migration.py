"""Library reprocessing preserves authored chemistry and migrates its users."""

import csv
import json
from pathlib import Path

import pytest
from rdkit import Chem

from pyPept.attachments import attachment_sites
from pyPept.library_quality import definition_hash
from pyPept.monomer_migration import (
    migrate_sequence,
    reprocess_library,
    reprocess_monomer,
    restored_sites,
)
from pyPept.monomer_store import _load_sdf, monomer_record


@pytest.mark.parametrize(
    "template,groups,slots,standalone",
    [
        (
            "[1*]N(C)[C@@H](CS[3*])C([2*])=O",
            {1: "[H]", 2: "[OH]", 3: "[H]"},
            {1: 1, 2: 2, 3: 4},
            "CN[C@@H](CS)C(=O)O",
        ),
        (
            "[1*]N([3*])[C@@H](CN1C(=O)C=C([4*])C1=O)C([2*])=O",
            {1: "[H]", 2: "[OH]", 3: "[H]", 4: "[H]"},
            {1: 1, 2: 2, 3: 3, 4: 4},
            "N[C@@H](CN1C(=O)C=CC1=O)C(=O)O",
        ),
        (
            "[1*]N([3*])/C(=C/C(=O)O)C(=O)O",
            {1: "[H]", 3: "[H]"},
            {1: 1, 3: 3},
            "N/C(=C/C(=O)O)C(=O)O",
        ),
    ],
)
def test_old_sites_keep_their_structure_and_roles(template, groups, slots, standalone):
    old = monomer_record(template, "Example", groups)
    result = reprocess_monomer(old)
    assert result.slot_mapping == slots
    assert Chem.MolToSmiles(restored_sites(result.molecule)[0]) == Chem.MolToSmiles(
        Chem.MolFromSmiles(standalone)
    )
    if "N1C(=O)" in template:
        assert len(attachment_sites(result.molecule)) == 4  # no second maleimide handle


def test_every_bundled_definition_is_stable_under_reprocessing_and_atom_order():
    for symbol, molecule in _load_sdf()[1].items():
        original = definition_hash(molecule)
        result = reprocess_monomer(molecule)
        assert definition_hash(result.molecule) == original, symbol
        assert all(old == new for old, new in result.slot_mapping.items()), symbol
        reversed_mol = Chem.RenumberAtoms(
            molecule, list(reversed(range(molecule.GetNumAtoms())))
        )
        for key in molecule.GetPropNames():
            reversed_mol.SetProp(key, molecule.GetProp(key))
        assert (
            definition_hash(reprocess_monomer(reversed_mol).molecule) == original
        ), symbol


def test_prochiral_attachment_keeps_the_authored_bonded_configuration():
    template = "[1*]N([3*])[C@@H]([2*])C(C)C"
    original = monomer_record(
        template,
        "ValAryl",
        {1: "[H]", 2: "[H]", 3: "[H]"},
        {1: "backbone_n", 2: "sp3_c_anchor", 3: "backbone_n_mod"},
    )
    processed = reprocess_monomer(original).molecule
    assert Chem.MolToSmiles(processed) == Chem.MolToSmiles(Chem.MolFromSmiles(template))
    # The standalone molecule has no stereocenter; that weaker check misses inversion.
    assert Chem.MolToSmiles(restored_sites(processed)[0]) == "CC(C)CN"


def test_bundled_csv_rebuild_retains_every_definition_and_subtype(tmp_path):
    from pyPept.interfaces import monomer_pipeline

    data = Path(monomer_pipeline.__file__).parents[1] / "data"
    rows, errors = monomer_pipeline.derive_monomers(data / "monomers.csv", rebuild=True)
    assert errors == []
    target = tmp_path / "monomers.sdf"
    count, errors = monomer_pipeline.write_sdf(rows, target)
    assert (count, errors) == (1128, [])
    expected = {name: definition_hash(mol) for name, mol in _load_sdf()[1].items()}
    actual = {
        mol.GetProp("symbol"): definition_hash(mol)
        for mol in Chem.SDMolSupplier(str(target), removeHs=False)
    }
    assert actual == expected


def test_migration_writes_complete_data_and_rebinds_nested_connections(tmp_path):
    old = tmp_path / "original.sdf"
    legacy = monomer_record(
        "[1*]N(C)[C@@H](CS[3*])C([2*])=O", "meC", {1: "[H]", 2: "[OH]", 3: "[H]"}
    )
    with Chem.SDWriter(str(old)) as writer:
        writer.write(legacy)
        writer.write(_load_sdf()[1]["C"])
    target = tmp_path / "migrated"
    manifest = reprocess_library(old, target)
    mappings = {name: record["slots"] for name, record in manifest["records"].items()}
    assert mappings["meC"] == {1: 1, 2: 2, 3: 4}
    assert (
        migrate_sequence("meC.[C(3,4)]", old, target / "monomers.sdf", mappings)
        == "meC.[C(4,4)]"
    )
    cross_chain = migrate_sequence(
        "meC.!s(3,4)%C.!s", old, target / "monomers.sdf", mappings
    )
    assert "(4,4)" in cross_chain
    assert (
        json.loads((target / "site-migration.json").read_text())["records"]["meC"][
            "slots"
        ]["3"]
        == 4
    )
    with (target / "monomers.csv").open(newline="") as stream:
        row = next(csv.DictReader(stream))
    assert row["r3_leaving"] == "" and row["r4_leaving"] == "[H]"
    assert row["chem_types"] == "1:backbone_n,2:backbone_c,4:thiol"
    with pytest.raises(ValueError, match="new output directory"):
        reprocess_library(old, target)
    with pytest.raises(ValueError):
        migrate_sequence("meC.[C(3,4)]", old, target / "monomers.sdf", {})


def test_explicit_csv_rebuild_processes_caps_as_well_as_amino_acids(tmp_path):
    from pyPept.interfaces.monomer_pipeline import derive_monomers

    path = tmp_path / "monomers.csv"
    path.write_text(
        "token,input,name,type,chuckles,r2_leaving,r3_leaving\n"
        "Mpa,,,cap,C[C@H](S[3*])C([2*])=O,[OH],[H]\n"
    )
    rows, errors = derive_monomers(path, rebuild=True)
    assert errors == []
    assert rows[0]["r3_leaving"] == ""
    assert rows[0]["r4_leaving"] == "[H]"
    assert rows[0]["activation_policy"] == "canonical-sites-v1"
