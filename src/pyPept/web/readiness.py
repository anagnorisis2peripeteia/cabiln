"""Validate selected application data once at startup, never on a probe."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from pyPept import monomer_store
from pyPept.library_quality import AUDIT_VERSION, quality_manifest


_ASSETS = (
    "index.html", "builder.js", "builder.css", "theme.css",
    "project.js", "document.js",
    "register.html", "register.js", "register.css", "examples.json",
)


def check_readiness(static_dir):
    path = monomer_store.library_path()
    binding = monomer_store.library_binding()
    assets = []
    for name in _ASSETS:
        file = static_dir / name
        stat = file.stat()
        if not stat.st_size:
            raise ValueError("A required browser asset is empty")
        assets.append((str(file), stat.st_size, stat.st_mtime_ns))
    manifest_file = (
        Path(monomer_store.__file__).resolve().parent / "data" / "library-quality.json"
    )
    stat = manifest_file.stat()
    if not stat.st_size:
        raise ValueError("The library quality manifest is empty")
    manifest = quality_manifest()
    if (
        not isinstance(manifest, dict) or manifest.get("schema_version") != 1
        or manifest.get("audit_version") != AUDIT_VERSION
        or not isinstance(manifest.get("records"), dict) or not manifest["records"]
    ):
        raise ValueError("The library quality manifest has an unsupported format")
    for record in manifest["records"].values():
        if (
            not isinstance(record, dict)
            or not isinstance(record.get("definition_hash"), str)
            or not isinstance(record.get("issues"), list)
            or not isinstance(record.get("activation"), dict)
            or not isinstance(record["activation"].get("status"), str)
        ):
            raise ValueError("A library quality record has an invalid format")
    _validate(str(path), tuple(sorted(binding.items())), tuple(assets))
    return {"ready": True, "library_binding": binding}


@lru_cache(maxsize=2)
def _validate(library, binding, assets):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    import yaml

    from pyPept.interfaces.reaction_library import REACTIONS, REACTION_INDEX

    with rdBase.BlockLogs():
        records = list(Chem.SDMolSupplier(library, removeHs=False))
        if not records or any(record is None for record in records):
            raise ValueError("Selected monomer library has invalid or no records")
        # This exercises the same aliases, attachment definitions and named
        # record checks used by request parsing, including external libraries.
        table = monomer_store.monomer_table(library)
        if len(table) != len(records):
            raise ValueError("Selected library has unresolved definitions")
        if not REACTIONS or not REACTION_INDEX:
            raise ValueError("Reaction rules are empty")
        for entry in REACTIONS.values():
            if not entry.get("steps"):
                raise ValueError("A reaction has no steps")
            for step in entry["steps"]:
                if AllChem.ReactionFromSmarts(step) is None:
                    raise ValueError("A reaction step is invalid")
    data = Path(monomer_store.__file__).resolve().parent / "data"
    caps = yaml.safe_load((data / "cap_reactions.yaml").read_text())
    if not isinstance(caps, dict) or not caps:
        raise ValueError("Cap definitions are empty or invalid")
    examples = Path(assets[-1][0])
    if not isinstance(json.loads(examples.read_text(encoding="utf-8")), list):
        raise ValueError("Browser examples must contain a list")
    if dict(binding) != monomer_store.library_binding():
        raise ValueError("Library changed during startup validation")
