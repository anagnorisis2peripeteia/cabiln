"""Bounded Hypothesis campaigns and replay evidence shared by contract owners."""

from collections import Counter, deque
from functools import lru_cache
import hashlib
from importlib.metadata import distributions
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from hypothesis import note, settings
from hypothesis.database import DirectoryBasedExampleDatabase

PROFILE = os.environ.get("CABILN_FUZZ_PROFILE", "ci")
if PROFILE not in {"ci", "deep"}:
    raise ValueError("CABILN_FUZZ_PROFILE must be ci or deep")
ARTIFACTS = os.environ.get("CABILN_FUZZ_ARTIFACTS")
_counts = Counter()
_recent = deque(maxlen=8)
_nodeid = None


def fuzz_settings(examples=50, **kwargs):
    """Deep runs multiply each owner's bounded example budget by ten."""
    options = dict(
        max_examples=examples * (10 if PROFILE == "deep" else 1),
        deadline=30_000,
        print_blob=True,
    )
    if ARTIFACTS:
        options["database"] = DirectoryBasedExampleDatabase(
            Path(ARTIFACTS) / "hypothesis"
        )
    options.update(kwargs)
    return settings(**options)


def record(feature, **metadata):
    """Record observations, including repeated examples during shrinking/replay."""
    observation = {"feature": feature, **metadata}
    _counts[feature] += 1
    _recent.append(observation)
    note(json.dumps(observation, sort_keys=True, default=repr))
    if ARTIFACTS:
        # A watchdog/native crash can bypass pytest's final report.
        path = Path(ARTIFACTS) / "active.json"
        temporary = path.with_suffix(".tmp")
        temporary.write_text(
            json.dumps(
                {"test": _nodeid, "recent_observations": list(_recent)}, default=repr
            )
            + "\n",
            encoding="utf-8",
        )
        temporary.replace(path)


def start(nodeid, invocation):
    global _nodeid
    _nodeid = nodeid
    _counts.clear()
    _recent.clear()
    if ARTIFACTS:
        destination = Path(ARTIFACTS)
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "context.json").write_text(
            json.dumps(
                {
                    **_context(),
                    "profile": PROFILE,
                    "pytest_arguments": list(invocation),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )


@lru_cache(maxsize=1)
def _context():
    root = Path(__file__).resolve().parents[1]
    paths = sorted((root / "src").rglob("*"))
    paths += sorted((root / "tests").glob("*.py"))
    hashes = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
        if path.is_file() and "__pycache__" not in path.parts
    }
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True
    )
    return {
        "commit": commit.stdout.strip() if commit.returncode == 0 else None,
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": dict(
            sorted(
                (package.metadata["Name"], package.version)
                for package in distributions()
            )
        ),
        "file_sha256": hashes,
    }


def finish(nodeid, report, invocation):
    """Persist pytest's falsifying example and reproduction blob on failure."""
    if not ARTIFACTS:
        return
    destination = Path(ARTIFACTS)
    destination.mkdir(parents=True, exist_ok=True)
    result = {
        "test": nodeid,
        "profile": PROFILE,
        "pytest_arguments": list(invocation),
        "outcome": report.outcome,
        "duration_seconds": report.duration,
        "observation_counts": dict(sorted(_counts.items())),
        "recent_observations": list(_recent),
        "failure": report.longreprtext if report.failed else None,
        "captured_output": report.sections if report.failed else [],
    }
    name = hashlib.sha256(nodeid.encode()).hexdigest()[:16]
    (destination / (name + ".json")).write_text(
        json.dumps(result, indent=2, default=repr) + "\n", encoding="utf-8"
    )
    (destination / "active.json").unlink(missing_ok=True)
