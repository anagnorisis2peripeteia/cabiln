"""Shared monomer persistence for the CLI, conversion and web application."""

from __future__ import annotations

import os
import threading
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path

_SDF_PATH = Path(__file__).resolve().parent / "data" / "monomers.sdf"
_sdf_cache = {"mols": None, "by_abbr": None, "version": None}
_sdf_lock = threading.RLock()
_table_cache = OrderedDict()


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


def library_version(path=None, *, include_aliases=True):
    """Identify the SDF and, by default, its optional synonym definitions.

    Raw structures and recognition patterns do not depend on CSV aliases. They
    use the same file identity stamp with ``include_aliases=False``.
    """
    path = library_path() if path is None else Path(path).expanduser().resolve()

    def stamp(file):
        stat = file.stat()
        return stat.st_mtime_ns, stat.st_size, stat.st_ctime_ns, stat.st_ino

    structure_version = str(path), stamp(path)
    if not include_aliases:
        return structure_version
    aliases = path.with_name("monomers.csv")
    try:
        alias_stamp = stamp(aliases)
    except FileNotFoundError:
        alias_stamp = None
    return (*structure_version, alias_stamp)


def monomer_table(path=None):
    """Return detached definitions from one stable library version.

    Parsing, alias resolution, and structural cap-family discovery belong to
    the library. Callers can add local synthetic records or mutate their own
    molecules without changing later consumers. Keep at most two file versions.
    """
    from rdkit import Chem

    path = library_path() if path is None else Path(path).expanduser().resolve()
    with _sdf_lock:
        for _ in range(3):
            version = library_version(path)
            if version not in _table_cache:
                table = _read_monomer_table(str(path))
                if library_version(path) != version:
                    continue
                _table_cache[version] = table
                while len(_table_cache) > 2:
                    _table_cache.popitem(last=False)
            _table_cache.move_to_end(version)
            original = _table_cache[version]
            result = original.copy(deep=True)
            # Pandas' deep copy retains mutable objects inside object columns.
            # Each Sequence gets separate molecules, leaving lists, and aliases.
            result["m_romol"] = [Chem.Mol(mol) for mol in original["m_romol"]]
            result["m_Rgroups"] = [list(value) for value in original["m_Rgroups"]]
            result.attrs = deepcopy(original.attrs)
            return result
    raise ValueError("Monomer library changed repeatedly while loading definitions")


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
    _, molecules, by_abbr = _load_sdf_snapshot()
    return molecules, by_abbr


def _load_sdf_snapshot():
    """Return raw definitions and their verified SDF version as one snapshot.

    Explicit hydrogens remain present. An external atomic replacement during a
    read causes a retry; data are never cached under a version from another file.
    """
    from rdkit import Chem

    with _sdf_lock:
        for _ in range(3):
            version = library_version(include_aliases=False)
            if _sdf_cache["mols"] is not None and _sdf_cache["version"] == version:
                return version, _sdf_cache["mols"], _sdf_cache["by_abbr"]
            supplier = Chem.SDMolSupplier(version[0], removeHs=False)
            molecules = [mol for mol in supplier if mol is not None]
            if library_version(include_aliases=False) != version:
                continue
            names = _index_monomer_names(
                (
                    mol.GetProp("symbol") if mol.HasProp("symbol") else None,
                    mol.GetProp("m_abbr") if mol.HasProp("m_abbr") else None,
                )
                for mol in molecules
            )
            by_abbr = {name: molecules[index] for name, index in names.items()}
            _sdf_cache.update(mols=molecules, by_abbr=by_abbr, version=version)
            return version, molecules, by_abbr
    raise ValueError("Monomer library changed repeatedly while loading structures")


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
            library = monomer_table(path)
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
        _sdf_cache["version"] = None
        _table_cache.clear()


def _read_monomer_table(path):
    """Normalize structures, terminal cap families, and aliases once per version."""
    from collections import defaultdict

    from rdkit import Chem
    from rdkit.Chem import PandasTools

    sdf_file = path
    df_group = PandasTools.LoadSDF(sdf_file)
    symbols = list(df_group["symbol"])
    abbreviations = list(df_group.get("m_abbr", [None] * len(symbols)))
    name_index = _index_monomer_names(zip(symbols, abbreviations))
    df_group["m_abbr"] = [
        abbr.strip() if isinstance(abbr, str) and abbr.strip() else symbol
        for symbol, abbr in zip(symbols, abbreviations)
    ]

    rg_col = "m_Rgroups"
    df_group[rg_col] = df_group[rg_col].astype(object)
    for idx in df_group.index:
        change = df_group[rg_col][idx].split(",")
        df_group.at[idx, rg_col] = [
            None if v.strip() in ("None", "") else v.strip() for v in change
        ]

    df_group = df_group.set_index("symbol")
    df_group = df_group.rename(columns={"ROMol": "m_romol"})

    # Propagate m_chem_types from SDF into each m_romol as an rdkit property
    # so downstream infer_chem_type can opt-in to SDF-declared chem types
    # (e.g. backbone_c_red for reduced-amide peptidomimetics) instead of
    # relying solely on structural heuristics.
    if "m_chem_types" in df_group.columns:
        for sym in df_group.index:
            mol = df_group.at[sym, "m_romol"]
            ct = df_group.at[sym, "m_chem_types"]
            if mol is not None and isinstance(ct, str) and ct:
                mol.SetProp("m_chem_types", ct)

    # --- Structural degeneracy detection for paired caps ---
    # Group cap monomers by their core structure (dummy atoms stripped).
    # Paired caps like Bn_/_Bn are structurally identical minus dummies;
    # we derive the base name and store it as a degenerate alias.
    caps = df_group[df_group["m_type"] == "cap"]
    core_map = {}  # canonical SMILES -> list of symbol names
    for sym in caps.index:
        mol = caps.at[sym, "m_romol"]
        if mol is None:
            continue
        # Strip dummy atoms (atomic num == 0) to get core
        emol = Chem.RWMol(mol)
        dummies = [a.GetIdx() for a in emol.GetAtoms() if a.GetAtomicNum() == 0]
        for idx in sorted(dummies, reverse=True):
            emol.RemoveAtom(idx)
        try:
            core_smi = Chem.MolToSmiles(emol.GetMol())
        except Exception:
            continue
        if core_smi not in core_map:
            core_map[core_smi] = []
        core_map[core_smi].append(sym)

    # For each group with >1 member, derive the base name and add alias
    degen_aliases = {}  # base_name -> list of variant symbols
    for smi, variants in core_map.items():
        if len(variants) < 2:
            continue
        # Derive base name: strip leading/trailing '_' from each variant
        bases = set()
        for v in variants:
            base = v.strip("_")
            bases.add(base)
        if len(bases) == 1:
            base_name = bases.pop()
            # Only create alias if base_name is not already an entry
            if base_name not in name_index:
                degen_aliases[base_name] = variants

    # Store aliases as a module-level accessible dict on the DataFrame
    df_group.attrs["_degen_aliases"] = degen_aliases

    # Stored symbols and abbreviations take precedence over CSV aliases.
    synonyms = {
        name: symbols[index]
        for name, index in name_index.items()
        if name != symbols[index]
    }
    csv_aliases = defaultdict(set)
    csv_path = os.path.join(os.path.dirname(path), "monomers.csv")
    if os.path.isfile(csv_path):
        import csv as _csv

        with open(csv_path, newline="", encoding="utf-8") as _f:
            for row in _csv.DictReader(_f):
                tok = row.get("token", "").strip()
                syns = row.get("synonyms", "")
                if tok and syns and tok in df_group.index:
                    for alias in syns.split(","):
                        alias = alias.strip()
                        if alias and alias not in name_index:
                            csv_aliases[alias].add(tok)
    ambiguous = {}
    for alias, targets in csv_aliases.items():
        if alias in degen_aliases:
            targets.update(degen_aliases[alias])
        if len(targets) == 1:
            synonyms[alias] = next(iter(targets))
        else:
            ambiguous[alias] = tuple(sorted(targets))
    df_group.attrs["_synonyms"] = synonyms
    df_group.attrs["_ambiguous_aliases"] = ambiguous

    return df_group
