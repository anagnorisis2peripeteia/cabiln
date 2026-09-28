#!/usr/bin/env python3
"""Bounded, offline benchmark of independently annotated molecular inputs.

Each case runs in a fresh process with CPU, wall-time and output-file limits.
Reference atom sets and bonds come from the fixture, never from the importer.
The optional external adapter is isolated and its missing evidence stays null.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _canonical(molecule):
    from rdkit import Chem

    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def _identity(source, product):
    from rdkit import Chem

    if product is None:
        return None, None
    exact = _canonical(source) == _canonical(product)
    source_connectivity, product_connectivity = Chem.Mol(source), Chem.Mol(product)
    Chem.RemoveStereochemistry(source_connectivity)
    Chem.RemoveStereochemistry(product_connectivity)
    # isomericSmiles=False also discards isotope labels. Stereo inference may
    # add an unspecified stereocentre, but must not change isotopes or charge.
    compatible = _canonical(source_connectivity) == _canonical(product_connectivity)
    compatible &= product.HasSubstructMatch(source, useChirality=True)
    return exact, bool(compatible)


def _validate(case):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        raise ValueError(f"Invalid reference structure: {case['id']}")
    for alternative in case["expected"]["alternatives"]:
        owners = {}
        for index, residue in enumerate(alternative["residues"]):
            for atom in residue["atoms"]:
                if atom in owners:
                    raise ValueError(f"Overlapping reference owners: {case['id']}")
                owners[atom] = index
            if not set(residue["required_slots"].values()) <= set(residue["atoms"]):
                raise ValueError(f"Reference attachment outside residue: {case['id']}")
        if set(owners) != set(range(molecule.GetNumAtoms())):
            raise ValueError(f"Incomplete reference ownership: {case['id']}")
        boundaries = {
            tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
            for bond in molecule.GetBonds()
            if owners[bond.GetBeginAtomIdx()] != owners[bond.GetEndAtomIdx()]
        }
        if boundaries != {tuple(sorted(edge)) for edge in alternative["connections"]}:
            raise ValueError(f"Reference boundaries disagree: {case['id']}")
    return molecule


def _evaluate(case, value, response):
    from rdkit import Chem

    reference = Chem.MolFromSmiles(case["smiles"])
    source = Chem.MolFromSmiles(value)
    product = None
    if response.get("rebuilt_smiles"):
        product = Chem.MolFromSmiles(response["rebuilt_smiles"])
    elif response.get("helm"):
        # This independent verifier supports only RDKit's HELM subset. Failure
        # to parse an external monomer is unavailable evidence, not proof of error.
        try:
            product = Chem.MolFromHELM(response["helm"])
        except (ValueError, RuntimeError):
            pass
    exact, compatible = _identity(source, product)
    metrics = {
        "exact_molecule": exact,
        "stereo_compatible": compatible,
        "atom_coverage": None,
        "expected_partition": None,
        "attachment_semantics": None,
        "known_retention": None,
        "emitted_ownership": response.get("emitted_ownership"),
        "emitted_recognition": response.get("emitted_recognition"),
        "accepted": None,
    }
    assignments = response.get("assignments")
    if assignments is None:
        return metrics
    metrics["atom_coverage"] = Counter(
        atom for assignment in assignments for atom in assignment["source_atoms"]
    ) == Counter(range(source.GetNumAtoms()))
    metrics["expected_partition"] = False
    mappings = source.GetSubstructMatches(
        reference, useChirality=True, uniquify=False, maxMatches=4096
    )
    if len(mappings) == 4096:
        raise ValueError("Reference symmetry exceeded the scoring limit")
    actual = {frozenset(a["source_atoms"]): a for a in assignments}
    actual_edges = {
        tuple(sorted(edge)) for edge in response.get("source_connections", [])
    }
    for alternative in case["expected"]["alternatives"]:
        for mapping in mappings:
            wanted = [
                frozenset(mapping[a] for a in r["atoms"])
                for r in alternative["residues"]
            ]
            if set(wanted) != set(actual) or len(wanted) != len(assignments):
                continue
            metrics["expected_partition"] = True
            known = all(
                residue["recognized"] is None
                or actual[atoms]["recognized"] is residue["recognized"]
                for atoms, residue in zip(wanted, alternative["residues"])
            )
            slots = all(
                dict(actual[atoms]["attachments"]).get(int(slot)) == mapping[atom]
                for atoms, residue in zip(wanted, alternative["residues"])
                for slot, atom in residue["required_slots"].items()
            )
            edges = {
                tuple(sorted(mapping[atom] for atom in edge))
                for edge in alternative["connections"]
            }
            metrics["known_retention"] = metrics["known_retention"] is True or known
            metrics["attachment_semantics"] = metrics[
                "attachment_semantics"
            ] is True or (slots and edges == actual_edges)
            if known and slots and edges == actual_edges:
                allowed_identity = exact or (
                    case["expected"]["allow_stereo_inference"]
                    and compatible
                    and response.get("inferred_stereo")
                )
                metrics["accepted"] = bool(
                    allowed_identity
                    and metrics["atom_coverage"]
                    and metrics["emitted_ownership"] is True
                    and metrics["emitted_recognition"] is True
                )
                return metrics
    metrics["accepted"] = False
    return metrics


def _emitted_witness(value, sequence, assembly, assignments):
    """Check the reported partition against the actual emitted occurrences."""
    from rdkit import Chem

    source = Chem.MolFromSmiles(value)
    product = assembly.get_molecule("ROMol")
    owners = assembly.get_residue_atom_map()
    synthetic = set(sequence._synthetic_aa_smiles) | set(sequence._synthetic_caps)
    indexes = [item["residue_index"] for item in assignments]
    if sorted(indexes) != list(range(len(sequence.s_monomers))):
        return {"emitted_ownership": False, "emitted_recognition": False}
    recognized = all(
        item["recognized"]
        == (sequence.s_monomers[item["residue_index"]]["m_abbr"] not in synthetic)
        for item in assignments
    )
    mappings = product.GetSubstructMatches(
        source, useChirality=True, uniquify=False, maxMatches=4096
    )
    consistent = len(mappings) < 4096 and any(
        all(
            set(owners.get(item["residue_index"], []))
            == {mapping[atom] for atom in item["source_atoms"]}
            for item in assignments
        )
        for mapping in mappings
    )
    return {"emitted_ownership": consistent, "emitted_recognition": recognized}


def _worker(request_path, result_path):
    request = json.loads(request_path.read_text())
    limits = request["case"]["limits"]
    cpu = limits["cpu_seconds"]
    resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu + 1))
    resource.setrlimit(
        resource.RLIMIT_FSIZE, (limits["output_bytes"], limits["output_bytes"])
    )
    from rdkit import Chem

    from pyPept.interfaces.cli_monomer import register_monomer
    from pyPept.molecule import Molecule
    from pyPept.monomer_store import library_path
    from pyPept.sequence import Sequence
    from pyPept.smiles import convert_smiles

    selected = library_path()
    if request["case"].get("library_overlay"):
        selected = request_path.parent / "monomers.sdf"
        original = library_path()
        shutil.copyfile(original, selected)
        companion = original.with_suffix(".csv")
        if companion.is_file():
            shutil.copyfile(companion, selected.with_suffix(".csv"))
        os.environ["CABILN_MONOMER_LIBRARY"] = str(selected)
        for monomer in request["case"]["library_overlay"]:
            register_monomer(monomer["smiles"], monomer["symbol"], sdf_path=selected)
    response = {
        "status": "error",
        "engine": "cabiln",
        "library_sha256": hashlib.sha256(selected.read_bytes()).hexdigest(),
    }
    started = time.perf_counter()
    try:
        result = convert_smiles(request["smiles"])
        elapsed = time.perf_counter() - started
        sequence = Sequence(result.cabiln, warning_sink=lambda _message: None)
        assembly = Molecule(sequence)
        product = assembly.get_molecule("ROMol")
        assignments = [asdict(a) for a in result.assignments]
        by_index = {a.residue_index: dict(a.attachments) for a in result.assignments}
        connections = []
        for left, _, right, _, left_slot, right_slot in sequence.s_bonds:
            if left in by_index and right in by_index:
                connections.append(
                    [by_index[left][left_slot], by_index[right][right_slot]]
                )
        response.update(
            status="success",
            elapsed_seconds=elapsed,
            cabiln=result.cabiln,
            recognition_status=result.recognition_status,
            search_complete=result.search_complete,
            inferred_stereo=result.inferred_stereo,
            warnings=result.warnings,
            assignments=assignments,
            source_connections=connections,
            rebuilt_smiles=Chem.MolToSmiles(product),
            **_emitted_witness(request["smiles"], sequence, assembly, assignments),
        )
    except Exception as error:
        response.update(
            error=f"{type(error).__name__}: {error}",
            elapsed_seconds=time.perf_counter() - started,
        )
    usage = resource.getrusage(resource.RUSAGE_SELF)
    response["cpu_seconds"] = usage.ru_utime + usage.ru_stime
    response["peak_rss_bytes"] = usage.ru_maxrss * (
        1 if sys.platform == "darwin" else 1024
    )
    result_path.write_text(json.dumps(response, indent=2) + "\n")


def _source_snapshot(source_root):
    paths = sorted(
        path
        for path in (source_root / "src/pyPept").rglob("*")
        if path.suffix in {".py", ".sdf", ".csv", ".yaml"}
    )
    return {
        str(p.relative_to(source_root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in paths
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, default=HERE / "decomposition_cases.json")
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--permutations", type=int, default=0)
    parser.add_argument("--engine", choices=["cabiln", "mol-to-helm"], default="cabiln")
    parser.add_argument("--external-checkout", type=Path)
    parser.add_argument("--external-library", type=Path)
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument("--revision-label", default=None)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker-request", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-result", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_request:
        _worker(args.worker_request, args.worker_result)
        return 0
    if args.permutations < 0 or args.permutations > 3:
        parser.error("--permutations must be between 0 and 3")
    if args.engine == "mol-to-helm" and not args.external_checkout:
        parser.error("--external-checkout is required for mol-to-helm")
    source_root = args.source_root.resolve()
    if not (source_root / "src/pyPept").is_dir():
        parser.error("--source-root must contain src/pyPept")
    from rdkit import Chem, rdBase

    cases = json.loads(args.cases.read_text())["cases"]
    if args.case:
        missing = set(args.case) - {case["id"] for case in cases}
        if missing:
            parser.error(f"Unknown cases: {sorted(missing)}")
        cases = [case for case in cases if case["id"] in args.case]
    output = args.output or Path(
        tempfile.mkdtemp(prefix="cabiln-independent-benchmark-")
    )
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        parser.error("--output must be an empty directory; recorded runs are preserved")
    snapshot = _source_snapshot(source_root)
    report = {
        "engine": args.engine,
        "git_revision": args.revision_label
        or subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_root": str(source_root),
        "python": sys.version,
        "rdkit": rdBase.rdkitVersion,
        "platform": platform.platform(),
        "fixture_sha256": hashlib.sha256(args.cases.read_bytes()).hexdigest(),
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": snapshot,
        "memory_limit": "not enforced; peak RSS reported for CABILN workers",
        "rows": [],
    }
    if args.engine == "mol-to-helm":
        report["external_revision"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=args.external_checkout, text=True
        ).strip()
        report["external_library"] = (
            {
                "path": str(args.external_library.resolve()),
                "sha256": hashlib.sha256(
                    args.external_library.read_bytes()
                ).hexdigest(),
            }
            if args.external_library
            else "upstream bundled HELMCoreLibrary.json"
        )
    for case in cases:
        reference = _validate(case)
        for variant in range(args.permutations + 1):
            value = case["smiles"]
            if variant:
                seed = (17, 41, 73)[variant - 1]
                order = list(range(reference.GetNumAtoms()))
                random.Random(seed).shuffle(order)
                renamed = Chem.RenumberAtoms(reference, order)
                value = Chem.MolToSmiles(renamed, canonical=False)
            identifier = f"{case['id']}-{variant}"
            folder = output / identifier
            folder.mkdir(exist_ok=True)
            request = {"id": identifier, "smiles": value, "case": case}
            request_path, result_path = folder / "input.json", folder / "result.json"
            request_path.write_text(json.dumps(request, indent=2) + "\n")
            command = [sys.executable, str(Path(__file__).resolve())]
            if args.engine == "cabiln":
                command += [
                    "--worker-request",
                    str(request_path),
                    "--worker-result",
                    str(result_path),
                ]
                data = None
            else:
                command = [
                    sys.executable,
                    str(HERE / "external_mol_to_helm.py"),
                    "--checkout",
                    str(args.external_checkout.resolve()),
                ]
                if args.external_library:
                    command += ["--library", str(args.external_library.resolve())]
                data = json.dumps(request) + "\n"

            def caps():
                cpu = case["limits"]["cpu_seconds"]
                resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu + 1))
                size = case["limits"]["output_bytes"]
                resource.setrlimit(resource.RLIMIT_FSIZE, (size, size))

            started = time.perf_counter()
            with (
                (folder / "stdout.log").open("w") as stdout,
                (folder / "stderr.log").open("w") as stderr,
            ):
                process = subprocess.Popen(
                    command,
                    cwd=source_root,
                    env={**os.environ, "PYTHONPATH": str(source_root / "src")},
                    stdin=subprocess.PIPE,
                    stdout=stdout,
                    stderr=stderr,
                    text=True,
                    start_new_session=True,
                    preexec_fn=caps,
                )
                try:
                    process.communicate(data, timeout=case["limits"]["wall_seconds"])
                    if args.engine == "mol-to-helm" and process.returncode == 0:
                        result_path.write_text((folder / "stdout.log").read_text())
                    response = (
                        json.loads(result_path.read_text())
                        if result_path.is_file()
                        else {
                            "status": (
                                "resource_limit"
                                if process.returncode < 0
                                else "worker_error"
                            ),
                            "returncode": process.returncode,
                        }
                    )
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    response = {"status": "timeout"}
            row = {
                "id": identifier,
                "case": case["id"],
                "variant": variant,
                "wall_seconds": time.perf_counter() - started,
                "response": response,
            }
            row["metrics"] = (
                _evaluate(case, value, response)
                if response.get("status") == "success"
                else {"accepted": False}
            )
            if args.engine != "cabiln" and case.get("library_overlay"):
                row["library_note"] = (
                    "CABILN-only overlay was not added to the external engine; libraries differ."
                )
            report["rows"].append(row)
            (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
            print(
                json.dumps(
                    {k: row[k] for k in ("id", "wall_seconds", "metrics")},
                    sort_keys=True,
                ),
                flush=True,
            )
    report["source_unchanged"] = snapshot == _source_snapshot(source_root)
    (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Evidence: {output / 'results.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
