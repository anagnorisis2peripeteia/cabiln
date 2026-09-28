"""Shared monomer persistence for the CLI, conversion and web application."""

from __future__ import annotations

import os
import threading
from pathlib import Path

_SDF_PATH = Path(__file__).resolve().parent / "data" / "monomers.sdf"
_sdf_cache = {"mols": None, "by_abbr": None, "mtime": 0}
_sdf_lock = threading.RLock()


def library_path():
    """Return the selected default SDF, optionally outside the installed package."""
    override = os.environ.get("CABILN_MONOMER_LIBRARY")
    if override:
        path = Path(override).expanduser().resolve()
        if not path.is_file():
            raise ValueError(
                f"CABILN_MONOMER_LIBRARY must name an existing SDF file: {path}"
            )
        return path
    return Path(_SDF_PATH).resolve()


def _index_monomer_names(identifiers):
    """Map canonical symbols and nonempty abbreviations to unique record positions."""
    names = {}
    for index, (symbol, abbreviation) in enumerate(identifiers):
        if not isinstance(symbol, str) or not symbol or symbol != symbol.strip():
            raise ValueError(
                "Each monomer needs a nonempty symbol without surrounding whitespace"
            )
        record_names = {symbol}
        if isinstance(abbreviation, str) and abbreviation.strip():
            record_names.add(abbreviation.strip())
        for name in record_names:
            if name in names:
                raise ValueError(
                    f"Ambiguous monomer name '{name}' belongs to multiple SDF records"
                )
            names[name] = index
    return names


def _load_sdf():
    """Return monomers and abbreviation lookup, reloading when path or file changes."""
    from rdkit import Chem

    with _sdf_lock:
        path = library_path()
        stat = path.stat()
        mtime = (str(path), stat.st_mtime_ns, stat.st_size)
        if _sdf_cache["mols"] is not None and _sdf_cache["mtime"] == mtime:
            return _sdf_cache["mols"], _sdf_cache["by_abbr"]
        suppl = Chem.SDMolSupplier(str(path), removeHs=False)
        mols = []
        for mol in suppl:
            if mol is None:
                continue
            mols.append(mol)
        names = _index_monomer_names(
            (
                mol.GetProp("symbol") if mol.HasProp("symbol") else None,
                mol.GetProp("m_abbr") if mol.HasProp("m_abbr") else None,
            )
            for mol in mols
        )
        by_abbr = {name: mols[index] for name, index in names.items()}
        _sdf_cache["mols"] = mols
        _sdf_cache["by_abbr"] = by_abbr
        _sdf_cache["mtime"] = mtime
        return mols, by_abbr


def register_molecule(mol, sdf_path=None):
    """Append a validated monomer atomically to the default or an explicit SDF.

    Missing and empty files are supported. Symbols, abbreviations and active
    parser aliases cannot be shadowed. The lock covers both the collision check
    and replacement, so concurrent CLI and web registrations share one writer.
    """
    import shutil
    import tempfile

    from filelock import FileLock
    from rdkit import Chem

    from pyPept.sequence import get_monomer_info
    from pyPept.structure import require_supported_stereo

    require_supported_stereo(mol)
    path = library_path() if sdf_path is None else Path(sdf_path).expanduser().resolve()
    names = set()
    for prop in ("m_abbr", "symbol"):
        value = mol.GetProp(prop) if mol.HasProp(prop) else ""
        if not value or value != value.strip():
            raise ValueError(
                f"Monomer {prop} must be nonempty without surrounding whitespace"
            )
        names.add(value)

    with FileLock(str(path) + ".lock", timeout=10), _sdf_lock:
        exists = path.exists()
        existing = []
        if exists and path.stat().st_size:
            existing = list(Chem.SDMolSupplier(str(path), removeHs=False))
            if any(entry is None for entry in existing):
                raise ValueError(f"Cannot register into an invalid monomer SDF: {path}")
        reserved = {
            entry.GetProp(prop).strip()
            for entry in existing
            for prop in ("m_abbr", "symbol")
            if entry.HasProp(prop)
        }
        # Reuse the parser's CSV synonyms and structurally derived cap aliases.
        if existing:
            library = get_monomer_info(str(path))
            reserved.update(library.attrs.get("_synonyms", {}))
            reserved.update(library.attrs.get("_degen_aliases", {}))
            reserved.update(library.attrs.get("_ambiguous_aliases", {}))
        collisions = names & reserved
        if collisions:
            name = sorted(collisions)[0]
            raise ValueError(
                f"Monomer name '{name}' already exists in library (or as an alias)"
            )
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=path.parent,
                prefix=".monomers-",
                suffix=".sdf",
                delete=False,
            ) as handle:
                temporary = Path(handle.name)
                if exists:
                    with path.open("rb") as original:
                        shutil.copyfileobj(original, handle)
                        if original.tell():
                            original.seek(-1, os.SEEK_END)
                            if original.read(1) != b"\n":
                                handle.write(b"\n")
                handle.write(Chem.SDWriter.GetText(mol, kekulize=False).encode("utf-8"))
                handle.flush()
                os.fsync(handle.fileno())
            if exists:
                temporary.chmod(path.stat().st_mode & 0o777)
            os.replace(temporary, path)
            _invalidate_sdf()
            return len(existing) + 1
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


def _invalidate_sdf():
    with _sdf_lock:
        _sdf_cache["mols"] = None
        _sdf_cache["by_abbr"] = None
        _sdf_cache["mtime"] = 0.0
