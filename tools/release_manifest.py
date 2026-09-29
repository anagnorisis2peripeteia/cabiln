"""Record installed chemistry bindings, runtime versions and the exact wheel."""

import argparse
import json
import os
import platform
from hashlib import sha256
from importlib.metadata import distributions
from pathlib import Path

from pyPept.canonical import canonical_convention
from pyPept.monomer_store import library_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    wheels = list(args.wheel_dir.glob("*.whl"))
    if len(wheels) != 1:
        parser.error("wheel directory must contain exactly one release wheel")
    manifest = {
        "release": os.environ.get("CABILN_RELEASE") or "local",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "wheel": {
            "name": wheels[0].name,
            "sha256": sha256(wheels[0].read_bytes()).hexdigest(),
        },
        "packages": dict(
            sorted((dist.metadata["Name"], dist.version) for dist in distributions())
        ),
        "library_binding": library_binding(),
        "canonical": canonical_convention(),
    }
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
