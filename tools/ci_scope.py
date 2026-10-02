"""Keep documentation-only changes out of the chemistry and release jobs."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import PurePosixPath


def is_documentation(path):
    return path in {
        "README.md",
        "CONTEXT.md",
        "CLAUDE.md",
        "tools/README.md",
        "tests/AGENTS.md",
        "tests/browser/README.md",
    } or (path.startswith("docs/") and PurePosixPath(path).suffix == ".md")


def changed_paths(event_name, event, *, cwd=None):
    """Return both sides of renames; unknown comparison bases require full CI."""
    if event_name == "pull_request":
        pull = event.get("pull_request", {})
        before, after = pull.get("base", {}).get("sha"), pull.get("head", {}).get("sha")
        comparison = "..."
    elif event_name == "push":
        before, after = event.get("before"), event.get("after")
        comparison = ".."
    else:
        return None
    if not all(
        isinstance(sha, str) and re.fullmatch(r"[0-9a-f]{40}", sha) and sha != "0" * 40
        for sha in (before, after)
    ):
        return None
    result = subprocess.run(
        [
            "git",
            "diff",
            "--no-renames",
            "--name-only",
            "-z",
            f"{before}{comparison}{after}",
            "--",
        ],
        cwd=cwd,
        capture_output=True,
    )
    if result.returncode:
        return None
    return [
        path.decode("utf-8", errors="surrogateescape")
        for path in result.stdout.split(b"\0")
        if path
    ]


def classify(paths):
    if paths is None:
        return {"application": True, "documentation": True}
    return {
        "application": any(not is_documentation(path) for path in paths),
        "documentation": any(is_documentation(path) for path in paths),
    }


if __name__ == "__main__":
    with open(os.environ["GITHUB_EVENT_PATH"], encoding="utf-8") as stream:
        event = json.load(stream)
    scope = classify(changed_paths(os.environ["GITHUB_EVENT_NAME"], event))
    print(json.dumps(scope))
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as stream:
        for key, value in scope.items():
            stream.write(f"{key}={str(value).lower()}\n")
