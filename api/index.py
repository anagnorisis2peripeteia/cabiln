import sys
import pathlib

# Make pyPept and tools importable from the repo root
_root = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_root / 'src'))
sys.path.insert(0, str(_root / 'tools'))

from live_renderer import app  # noqa: F401 — Vercel discovers `app`
