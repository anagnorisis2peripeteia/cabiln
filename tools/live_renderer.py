#!/usr/bin/env python3
"""Compatibility launcher. Implementation lives in pyPept.web and pyPept.smiles."""
from importlib import import_module


def __getattr__(name):
    # Old scripts imported conversion helpers from this development launcher.
    # Keep those imports working without making the chemistry library import FastAPI.
    modules = ("pyPept.smiles", "pyPept.web.notation", "pyPept.web.drawing",
               "pyPept.web.schemas", "pyPept.web.app", "pyPept.web.rendering",
               "pyPept.web.monomers", "pyPept.web.builder", "pyPept.web.conversion")
    for module_name in modules:
        module = import_module(module_name)
        if hasattr(module, name):
            return getattr(module, name)
    raise AttributeError(f"{__name__!r} has no attribute {name!r}")


if __name__ == "__main__":
    from pyPept.web.app import main
    main()
