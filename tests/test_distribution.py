"""Exercise the release artifact without editable installs or checkout imports."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _run(*command, cwd):
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout
    return result.stdout


@pytest.mark.distribution
def test_built_wheel_runs_core_and_web_outside_checkout(tmp_path):
    """Build through the sdist, install declared deps, then use both public layers."""
    distribution_dir = tmp_path / "dist"
    _run(
        sys.executable,
        "-m",
        "build",
        "--outdir",
        str(distribution_dir),
        str(PROJECT_ROOT),
        cwd=tmp_path,
    )
    wheel = next(distribution_dir.glob("*.whl"))
    environment = tmp_path / "installed"
    _run(sys.executable, "-m", "venv", str(environment), cwd=tmp_path)
    scripts = environment / ("Scripts" if os.name == "nt" else "bin")
    python = str(scripts / ("python.exe" if os.name == "nt" else "python"))
    _run(
        python,
        "-I",
        "-m",
        "pip",
        "install",
        "--disable-pip-version-check",
        str(wheel),
        cwd=tmp_path,
    )

    _run(
        python,
        "-I",
        "-c",
        textwrap.dedent("""
            import importlib.util
            import sys
            from importlib.metadata import metadata
            from importlib.resources import files
            from pathlib import Path

            import pyPept
            from pyPept.interfaces.reaction_library import REACTIONS
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence
            from pyPept.smiles import convert_smiles, smiles_to_cabiln_core
            from rdkit import Chem

            assert Path(pyPept.__file__).is_relative_to(Path(sys.prefix))
            assert importlib.util.find_spec("fastapi") is None
            assert metadata("pyPept")["Description-Content-Type"] == "text/markdown"
            for name in (
                "monomers.sdf", "monomers.csv", "reactions.yaml",
                "cap_reactions.yaml", "matrix.txt", "total_SS.txt",
                "library-quality.json",
            ):
                assert files("pyPept.data").joinpath(name).read_bytes(), name
            assert REACTIONS
            assert callable(smiles_to_cabiln_core)
            expected = Chem.MolToSmiles(
                Chem.MolFromSmiles("N[C@@H](C)C(=O)NCC(=O)O")
            )
            for notation in ("A-G", "Ala-Gly"):
                molecule = Molecule(Sequence(notation)).get_molecule(fmt="ROMol")
                assert Chem.MolToSmiles(molecule) == expected
            conversion = convert_smiles(expected)
            assert conversion.recognition_status == "complete"
            assert [item.symbol for item in conversion.assignments] == ["A", "G"]
            rebuilt = Molecule(Sequence(conversion.cabiln)).get_molecule(fmt="ROMol")
            assert Chem.MolToSmiles(rebuilt) == expected
            """),
        cwd=tmp_path,
    )
    for entry_point in ("run_pyPept", "pyPept-BILN-validate", "pyPept-monomer-add"):
        executable = scripts / (entry_point + (".exe" if os.name == "nt" else ""))
        assert "usage:" in _run(str(executable), "--help", cwd=tmp_path).lower()

    _run(
        python,
        "-I",
        "-m",
        "pip",
        "install",
        "--disable-pip-version-check",
        f"{wheel}[web]",
        "httpx",
        cwd=tmp_path,
    )
    _run(
        python,
        "-I",
        "-c",
        textwrap.dedent("""
            from fastapi.testclient import TestClient
            from pyPept.web.app import app, create_app

            with TestClient(app) as client:
                response = client.get("/")
                assert response.status_code == 200
                assert "/static/builder.js" in response.text
                assert client.get("/capabilities").json() == {"registration": False}
                assert client.get("/register").status_code == 403
                for name in (
                    "theme.css", "builder.js", "builder.css", "register.js",
                    "register.css", "examples.json", "project.js", "document.js",
                ):
                    response = client.get("/static/" + name)
                    assert response.status_code == 200, name
                    assert response.content, name
                assert client.get("/health").json() == {"status": "ok"}
                assert client.get("/ready").status_code == 200
                assert client.get("/examples").json()
                response = client.post("/render", json={"cabiln": "A-G"})
                assert response.status_code == 200, response.text
                assert "<svg" in response.json()["svg"]
            with TestClient(create_app(allow_registration=True)) as client:
                assert client.get("/capabilities").json() == {"registration": True}
                assert client.get("/register").status_code == 200
            """),
        cwd=tmp_path,
    )
    executable = scripts / ("cabiln.exe" if os.name == "nt" else "cabiln")
    assert "usage:" in _run(str(executable), "--help", cwd=tmp_path).lower()
    _run(python, "-I", "-m", "pip", "check", cwd=tmp_path)
