"""Backups preserve custom structures and aliases, including failed restores."""

import json
import zipfile

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import library_snapshots as snapshots
from pyPept import monomer_store
from pyPept.web.app import create_app


@pytest.fixture
def library(tmp_path, monkeypatch):
    _, records = monomer_store._load_sdf()
    path = tmp_path / "monomers.sdf"
    path.write_text(Chem.SDWriter.GetText(records["G"], kekulize=False))
    monkeypatch.setattr(monomer_store, "_SDF_PATH", path)
    monomer_store._invalidate_sdf()
    yield path
    monomer_store._invalidate_sdf()


def test_snapshot_restore_keeps_both_files_and_pre_restore_copy(library):
    aliases = library.with_name("monomers.csv")
    aliases.write_text("symbol,synonyms\nG,my-glycine\n")
    before = {path.name: path.read_bytes() for path in (library, aliases)}
    directory = library.parent / "backups"
    saved = snapshots.snapshot_library(library, directory)
    assert saved.stat().st_mode & 0o777 == 0o600
    aliases.write_text("symbol,synonyms\nG,changed\n")
    previous_content = aliases.read_bytes()
    previous = snapshots.restore_library(saved, library, directory, offline=True)
    assert {p.name: p.read_bytes() for p in (library, aliases)} == before
    assert snapshots._read_snapshot(previous)["monomers.csv"] == previous_content
    assert saved.exists()


def test_restore_requires_stopped_readers(library):
    with pytest.raises(ValueError, match="Stop the service"):
        snapshots.restore_library("unused", library, library.parent / "backups")


def test_restore_removes_aliases_absent_in_original(library):
    directory = library.parent / "backups"
    saved = snapshots.snapshot_library(library, directory)
    aliases = library.with_name("monomers.csv")
    aliases.write_text("symbol,synonyms\nG,temporary\n")
    snapshots.restore_library(saved, library, directory, offline=True)
    assert not aliases.exists()


@pytest.mark.parametrize("failure", ["hash", "rules", "extra-entry"])
def test_corrupt_snapshot_rejected_without_changing_library(library, failure):
    directory = library.parent / "backups"
    original = library.read_bytes()
    saved = snapshots.snapshot_library(library, directory)
    with zipfile.ZipFile(saved) as archive:
        entries = {name: archive.read(name) for name in archive.namelist()}
    if failure == "hash":
        entries["monomers.sdf"] += b"corruption"
    elif failure == "rules":
        manifest = json.loads(entries["snapshot.json"])
        manifest["rules"]["reactions.yaml"] = "changed"
        entries["snapshot.json"] = json.dumps(manifest)
    else:
        entries["../escape"] = b"not allowed"
    corrupted = library.parent / "corrupted.zip"
    with zipfile.ZipFile(corrupted, "w") as archive:
        for name, content in entries.items():
            archive.writestr(name, content)
    with pytest.raises(ValueError):
        snapshots.restore_library(corrupted, library, directory, offline=True)
    assert library.read_bytes() == original
    assert len(list(directory.glob("*.zip"))) == 1


def test_partial_restore_failure_rolls_back_pair(library, monkeypatch):
    directory = library.parent / "backups"
    saved = snapshots.snapshot_library(library, directory)
    aliases = library.with_name("monomers.csv")
    aliases.write_text("symbol,synonyms\nG,current\n")
    original = library.read_bytes(), aliases.read_bytes()
    replace = snapshots._replace_file
    calls = 0

    def fail_second(path, content):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("Simulated alias write failure")
        return replace(path, content)

    monkeypatch.setattr(snapshots, "_replace_file", fail_second)
    with pytest.raises(OSError):
        snapshots.restore_library(saved, library, directory, offline=True)
    assert (library.read_bytes(), aliases.read_bytes()) == original


@pytest.mark.parametrize("backup_fails", [False, True])
def test_authenticated_registration_backs_up_before_mutation(
    library, monkeypatch, backup_fails
):
    directory = library.parent / "backups"
    original = library.read_bytes()
    monkeypatch.setenv("CABILN_LIBRARY_BACKUP_DIR", str(directory))
    token = "test-only-administrative-token-123456789"
    monkeypatch.setenv("CABILN_REGISTRATION_TOKEN", token)
    if backup_fails:

        def fail_backup(*args):
            raise OSError("Simulated unavailable backup storage")

        monkeypatch.setattr(snapshots, "_snapshot_locked", fail_backup)
    with TestClient(create_app(allow_registration=True)) as client:
        preview = client.post(
            "/preview_monomer",
            json={
                "smiles": "N[C@@H](C)C(=O)O",
            },
        ).json()
        payload = {key: preview[key] for key in ("chuckles", "chem_types", "leaving")}
        payload.update(abbr="SavedA", name="Saved alanine")
        response = client.post("/register_monomer", json=payload, auth=("admin", token))
        assert response.status_code == (500 if backup_fails else 200), response.text
    if backup_fails:
        assert library.read_bytes() == original
    else:
        (saved,) = directory.glob("*.zip")
        assert snapshots._read_snapshot(saved)["monomers.sdf"] == original
        assert len(list(Chem.SDMolSupplier(str(library)))) == 2
