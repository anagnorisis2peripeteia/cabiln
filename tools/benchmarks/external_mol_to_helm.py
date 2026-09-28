#!/usr/bin/env python3
"""Isolated JSON-lines adapter for an explicitly supplied mol-to-helm checkout.

This does not install or vendor Datagrok's code. Its public API returns HELM,
without a reconstruction or a supported source-atom ownership result. Therefore
``status=success`` records only upstream acceptance; ``verified`` stays false.
Run it using an interpreter with the checkout's requirements already available.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", required=True, type=Path)
    parser.add_argument("--library", type=Path, help="Optional HELM monomer JSON")
    args = parser.parse_args()
    logics = args.checkout.resolve() / "logics"
    if not (logics / "pipeline.py").is_file():
        parser.error("--checkout must contain logics/pipeline.py")
    sys.path.insert(0, str(logics))
    with redirect_stdout(sys.stderr):
        from pipeline import convert_molecules_batch

    library_json = args.library.read_text() if args.library else None
    for line in sys.stdin:
        if not line.strip():
            continue
        started = time.perf_counter()
        result = {
            "id": None,
            "engine": "datagrok-mol-to-helm",
            "verified": False,
            "ownership_status": "unavailable_from_public_api",
            "reconstruction_status": "not_assessed",
        }
        try:
            request = json.loads(line)
            result["id"] = request["id"]
            smiles = request["smiles"]
            if not isinstance(smiles, str):
                raise ValueError("smiles must be a string")
            with redirect_stdout(sys.stderr):
                success, output = convert_molecules_batch(
                    [smiles], library_json=library_json, input_type="smiles"
                )[0]
            result["status"] = "success" if success else "failed"
            result["helm" if success else "error"] = output
        except Exception as error:
            result["status"] = "error"
            result["error"] = f"{type(error).__name__}: {error}"
        result["elapsed_seconds"] = time.perf_counter() - started
        print(json.dumps(result, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
