#!/usr/bin/env python3
"""Execute the complete Python examples published in the README."""

from check_docs import ROOT, python_blocks

if __name__ == "__main__":
    examples = python_blocks((ROOT / "README.md").read_text(encoding="utf-8"))
    if not examples:
        raise SystemExit("README has no Python examples to check")
    for index, code in enumerate(examples, 1):
        exec(
            compile(code, f"README.md example {index}", "exec"),
            {"__name__": "__main__"},
        )
    print(f"Executed {len(examples)} README Python examples.")
