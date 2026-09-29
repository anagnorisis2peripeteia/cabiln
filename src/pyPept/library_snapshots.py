"""Versioned backups and offline restoration of an administrative SDF library."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
import uuid
import zipfile
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path

from filelock import FileLock

MAX_SNAPSHOT_BYTES = 128 * 1024 * 1024


def _rules():
    data = Path(__file__).resolve().parent / "data"
    return {
        name: sha256((data / name).read_bytes()).hexdigest()
        for name in ("reactions.yaml", "cap_reactions.yaml")
    }


def _snapshot_locked(library, directory):
    """Caller holds the same library writer lock as registration."""
    library, directory = Path(library), Path(directory)
    if not library.is_file():
        raise ValueError("Snapshot requires an existing monomer library")
    aliases = library.with_name("monomers.csv")
    files = {"monomers.sdf": library}
    if aliases.exists():
        files["monomers.csv"] = aliases
    if sum(path.stat().st_size for path in files.values()) > MAX_SNAPSHOT_BYTES:
        raise ValueError("Library exceeds the snapshot size limit")
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    manifest = {
        "format": "cabiln-library-snapshot",
        "version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "files": {},
        "rules": _rules(),
    }
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = directory / f"library-{stamp}-{uuid.uuid4().hex[:8]}.zip"
    with tempfile.NamedTemporaryFile(dir=directory, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as archive:
            for name, path in files.items():
                content = path.read_bytes()
                manifest["files"][name] = sha256(content).hexdigest()
                archive.writestr(name, content)
            archive.writestr("snapshot.json", json.dumps(manifest, sort_keys=True))
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        return destination
    finally:
        temporary.unlink(missing_ok=True)


def snapshot_library(library, directory):
    """Capture the SDF and companion aliases under the registration writer lock."""
    from pyPept.monomer_store import _sdf_lock

    library = Path(library).expanduser().resolve()
    with FileLock(str(library) + ".lock", timeout=10), _sdf_lock:
        return _snapshot_locked(library, Path(directory).expanduser().resolve())


def _read_snapshot(snapshot):
    allowed = {"snapshot.json", "monomers.sdf", "monomers.csv"}
    with zipfile.ZipFile(snapshot) as archive:
        entries = archive.infolist()
        names = [entry.filename for entry in entries]
        if (
            len(names) != len(set(names))
            or not set(names).issubset(allowed)
            or not {"snapshot.json", "monomers.sdf"}.issubset(names)
            or sum(entry.file_size for entry in entries) > MAX_SNAPSHOT_BYTES
        ):
            raise ValueError("Invalid or oversized library snapshot")
        manifest = json.loads(archive.read("snapshot.json"))
        if (
            manifest.get("format") != "cabiln-library-snapshot"
            or manifest.get("version") != 1
            or manifest.get("rules") != _rules()
        ):
            raise ValueError(
                "Snapshot format or reaction rules differ from this release"
            )
        files = {name: archive.read(name) for name in names if name != "snapshot.json"}
        if manifest.get("files") != {
            name: sha256(content).hexdigest() for name, content in files.items()
        }:
            raise ValueError("Snapshot contents do not match their recorded hashes")
        return files


def restore_library(snapshot, library, backup_directory, *, offline=False):
    """Restore while service readers are stopped; retain the pre-restore snapshot.

    The SDF and CSV are individually replaced atomically. They are a pair, so
    live readers must be stopped for the entire operation, including rollback.
    """
    from rdkit import Chem

    from pyPept.monomer_store import _invalidate_sdf, _read_monomer_table, _sdf_lock

    if not offline:
        raise ValueError("Stop the service and explicitly select offline restoration")
    files = _read_snapshot(snapshot)
    library = Path(library).expanduser().resolve()
    aliases = library.with_name("monomers.csv")
    library.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=library.parent) as temporary:
        staged = Path(temporary)
        for name, content in files.items():
            (staged / name).write_bytes(content)
        molecules = list(
            Chem.SDMolSupplier(str(staged / "monomers.sdf"), removeHs=False)
        )
        if not molecules or any(mol is None for mol in molecules):
            raise ValueError("Snapshot contains an invalid or empty SDF")
        _read_monomer_table(str(staged / "monomers.sdf"))
        with FileLock(str(library) + ".lock", timeout=10), _sdf_lock:
            previous = (
                _snapshot_locked(library, backup_directory)
                if library.exists()
                else None
            )
            originals = {
                library: library.read_bytes() if library.exists() else None,
                aliases: aliases.read_bytes() if aliases.exists() else None,
            }
            try:
                _replace_file(library, files["monomers.sdf"])
                _replace_file(aliases, files.get("monomers.csv"))
            except BaseException:
                for path, content in originals.items():
                    _replace_file(path, content)
                raise
            finally:
                _invalidate_sdf()
            return previous


def _replace_file(path, content):
    if content is None:
        path.unlink(missing_ok=True)
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        if path.exists():
            temporary.chmod(path.stat().st_mode & 0o777)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("snapshot", "restore"))
    parser.add_argument("--library", required=True)
    parser.add_argument("--directory", required=True, help="Backup directory")
    parser.add_argument("--snapshot", help="Snapshot ZIP to restore")
    parser.add_argument("--offline", action="store_true", help="Service is stopped")
    args = parser.parse_args()
    if args.action == "snapshot":
        print(snapshot_library(args.library, args.directory))
    else:
        if not args.snapshot or not args.offline:
            parser.error("restore requires --snapshot and --offline")
        previous = restore_library(
            args.snapshot, args.library, args.directory, offline=True
        )
        print(f"Library restored; previous snapshot: {previous or 'new library'}")


if __name__ == "__main__":
    main()
