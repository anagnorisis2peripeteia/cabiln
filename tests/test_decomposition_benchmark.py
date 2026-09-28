"""Negative controls for the independent benchmark's acceptance contract.

A molecular round-trip alone, a plausible partition alone, or correct residue
names with wrong attachment atoms must each fail the benchmark independently.
These controls do not call the converter to manufacture their expected answers.
"""

import importlib.util
import json
from pathlib import Path

import pytest

BENCHMARKS = Path(__file__).resolve().parents[1] / "tools" / "benchmarks"
SPEC = importlib.util.spec_from_file_location(
    "decomposition_benchmark", BENCHMARKS / "decomposition.py"
)
BENCHMARK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCHMARK)
CASES = json.loads((BENCHMARKS / "decomposition_cases.json").read_text())["cases"]
ALA_GLY = next(case for case in CASES if case["id"] == "ala_gly")


def _independent_answer():
    return {
        "rebuilt_smiles": "N[C@@H](C)C(=O)NCC(=O)O",
        "assignments": [
            {
                "source_atoms": [0, 1, 2, 3, 4],
                "attachments": [[1, 0], [2, 3]],
                "recognized": True,
                "residue_index": 0,
            },
            {
                "source_atoms": [5, 6, 7, 8, 9],
                "attachments": [[1, 5], [2, 7]],
                "recognized": True,
                "residue_index": 1,
            },
        ],
        "source_connections": [[3, 5]],
        "emitted_ownership": True,
        "emitted_recognition": True,
    }


def test_reference_annotations_cover_every_atom_and_only_actual_boundary_bonds():
    for case in CASES:
        BENCHMARK._validate(case)


@pytest.mark.parametrize("failure", [None, "opaque", "reversed", "slots"])
def test_molecule_partition_and_slots_are_separate_requirements(failure):
    response = _independent_answer()
    if failure == "opaque":
        response["assignments"] = [
            {
                "source_atoms": list(range(10)),
                "attachments": [],
                "recognized": False,
            }
        ]
        response["source_connections"] = []
    elif failure == "reversed":
        response["rebuilt_smiles"] = "NCC(=O)N[C@@H](C)C(=O)O"
    elif failure == "slots":
        response["assignments"][0]["attachments"] = [[1, 3], [2, 0]]
    scored = BENCHMARK._evaluate(ALA_GLY, ALA_GLY["smiles"], response)
    assert scored["accepted"] is (failure is None)
    assert scored["atom_coverage"] is True
    assert scored["exact_molecule"] is (failure != "reversed")
    assert scored["expected_partition"] is (failure != "opaque")
    if failure == "slots":
        assert scored["attachment_semantics"] is False


def test_external_success_without_a_witness_does_not_pass_partition_acceptance():
    response = {"status": "success", "helm": "PEPTIDE1{A.G}$$$$V2.0"}
    scored = BENCHMARK._evaluate(ALA_GLY, ALA_GLY["smiles"], response)
    assert scored["exact_molecule"] is True
    assert scored["accepted"] is None
    assert scored["atom_coverage"] is None
    assert scored["expected_partition"] is None


def test_allowed_stereo_inference_cannot_hide_a_changed_isotope():
    case = next(case for case in CASES if case["id"] == "unspecified_ala_gly")
    response = _independent_answer()
    response["inferred_stereo"] = True
    response["rebuilt_smiles"] = "N[C@@H]([13CH3])C(=O)NCC(=O)O"
    scored = BENCHMARK._evaluate(case, case["smiles"], response)
    assert scored["expected_partition"] is True
    assert scored["stereo_compatible"] is False
    assert scored["accepted"] is False


def test_correct_atom_sets_with_wrong_emitted_occurrence_indices_fail():
    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence

    sequence = Sequence("A-G")
    assembly = Molecule(sequence)
    assignments = _independent_answer()["assignments"]
    witnessed = BENCHMARK._emitted_witness(
        ALA_GLY["smiles"], sequence, assembly, assignments
    )
    assert witnessed["emitted_ownership"] is True
    assignments[0]["residue_index"], assignments[1]["residue_index"] = 1, 0
    witnessed = BENCHMARK._emitted_witness(
        ALA_GLY["smiles"], sequence, assembly, assignments
    )
    assert witnessed["emitted_ownership"] is False
