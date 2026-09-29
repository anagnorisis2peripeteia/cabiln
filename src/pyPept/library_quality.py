"""Versioned software-compatibility facts for bundled monomer definitions.

The audit records activation discrepancies, metadata differences and unspecified
stereo; it does not certify scientific identity or reaction feasibility. Runtime
lookups only fingerprint a definition and read the packaged baseline. Activation
is performed by the explicit audit command and the library regression test.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
from importlib.resources import files
from pathlib import Path

from rdkit import Chem, rdBase

from pyPept.leaving_groups import restore_leaving_groups
from pyPept.structure import require_supported_stereo

AUDIT_VERSION = "cabiln-library-quality-v1"


def _property(molecule, name):
    return molecule.GetProp(name).strip() if molecule.HasProp(name) else ""


def _leaving(molecule):
    return {
        slot: value.strip()
        for slot, value in enumerate(_property(molecule, "m_Rgroups").split(","), 1)
        if value.strip().lower() not in ("", "none")
    }


def _digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return sha256(encoded.encode("utf-8")).hexdigest()


def definition_hash(molecule):
    """Fingerprint structure, named identity, numbered sites and their metadata.

    Coordinates and atom traversal do not identify a definition. Unknown stereo
    groups cannot accidentally acquire the fingerprint of an absolute template.
    This convention is tested on the supported RDKit versions; a changed hash
    safely makes the definition unreviewed until its baseline is checked again.
    """
    require_supported_stereo(molecule)
    return _digest(
        {
            "template": Chem.MolToSmiles(molecule, isomericSmiles=True),
            "identity": {
                key: _property(molecule, key)
                for key in ("symbol", "m_abbr", "m_type", "m_subtype")
            },
            "leaving": _leaving(molecule),
            "declared_chemistry": sorted(
                part.strip()
                for part in _property(molecule, "m_chem_types").split(",")
                if part.strip()
            ),
        }
    )


@lru_cache(maxsize=1)
def _manifest():
    resource = files("pyPept.data").joinpath("library-quality.json")
    return json.loads(resource.read_text(encoding="utf-8"))


def quality_manifest():
    """Return a detached copy of the packaged review baseline."""
    return deepcopy(_manifest())


def quality_for_monomer(molecule):
    """Return JSON-safe facts only for an exactly matching reviewed definition.

    Newly registered/custom/changed records remain available to every consumer,
    with ``unreviewed`` status. No activation work or library mutation happens
    here. ``no_known_exception`` is not a claim of scientific validation.
    """
    symbol = _property(molecule, "symbol") or _property(molecule, "m_abbr")
    try:
        fingerprint = definition_hash(molecule)
    except (ValueError, RuntimeError):
        fingerprint = None
    record = _manifest()["records"].get(symbol)
    matched = record is not None and record["definition_hash"] == fingerprint
    issues = deepcopy(record["issues"]) if matched else []
    return {
        "status": (
            ("review_required" if issues else "no_known_exception")
            if matched
            else "unreviewed"
        ),
        "issues": issues,
        "definition_hash": fingerprint,
        "audit_version": AUDIT_VERSION,
        "activation_status": record["activation"]["status"] if matched else None,
    }


def _issue(code, message, **facts):
    return {"code": code, "message": message, **facts}


def audit_monomer(molecule):
    """Measure one definition without altering its stored structure or slots."""
    from pyPept.attachments import attachment_sites
    from pyPept.interfaces.monomer_pipeline import pre_activate

    mol = Chem.Mol(molecule)
    leaving = _leaving(mol)
    slots = {
        atom.GetIsotope(): atom.GetNeighbors()[0].GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetDegree() == 1
    }
    legacy = all(slot in slots for slot in (1, 2, 3)) and slots[3] != slots[1]
    issues = []
    if legacy:
        issues.append(
            _issue(
                "legacy_r3_sidechain",
                "Stored R3 is a sidechain site; existing slot numbers are preserved.",
                slots=sorted(slot for slot in slots if slot >= 3),
            )
        )
    missing = sorted(set(slots) - set(leaving))
    if missing:
        issues.append(
            _issue(
                "missing_leaving_metadata",
                "Some sites lack leaving metadata; current restoration assumes H.",
                slots=missing,
            )
        )
    groups = [leaving.get(slot) for slot in range(1, max(slots, default=0) + 1)]
    differences = [
        {
            "slot": site["slot"],
            "declared": site["declared_chem_type"],
            "effective": site["chem_type"],
        }
        for site in attachment_sites(mol, groups)
        if site["declared_chem_type"]
        and site["declared_chem_type"] != site["chem_type"]
    ]
    if differences:
        issues.append(
            _issue(
                "chemistry_declaration_difference",
                "Declared and effective site chemistry differ; "
                "this is not proof of error.",
                sites=differences,
            )
        )
    expected = Chem.Mol(mol)
    expected_leaving = dict(leaving)
    if legacy:
        # Comparison only: never change the library or user-visible numbering.
        for atom in expected.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() >= 3:
                atom.SetIsotope(atom.GetIsotope() + 1)
        expected_leaving = {
            slot + 1 if slot >= 3 else slot: value for slot, value in leaving.items()
        }
    expected_template = Chem.MolToSmiles(expected, isomericSmiles=True)
    activation = {
        "status": "restoration_failed",
        "expected_template": expected_template,
        "expected_leaving": {str(k): v for k, v in sorted(expected_leaving.items())},
    }
    try:
        restored = restore_leaving_groups(mol, leaving)
        unassigned = sum(
            label == "?"
            for _, label in Chem.FindMolChiralCenters(
                restored,
                includeUnassigned=True,
                useLegacyImplementation=False,
            )
        )
        if unassigned:
            issues.append(
                _issue(
                    "unspecified_tetrahedral_stereo",
                    "The stored standalone structure leaves "
                    "tetrahedral stereo unspecified.",
                    count=unassigned,
                )
            )
        activation["status"] = "activation_failed"
        activated = pre_activate(Chem.MolToSmiles(restored, isomericSmiles=True))
        actual = Chem.MolFromSmiles(activated.chuckles)
        actual_template = Chem.MolToSmiles(actual, isomericSmiles=True)
        changes = []
        if actual_template != expected_template:
            changes.append("template")
        if activated.leaving != expected_leaving:
            changes.append("leaving_groups")
        if activated.chem_types.get(3) not in (None, "backbone_n_mod"):
            changes.append("reserved_r3_chemistry")
        activation.update(
            status=(
                "mismatch" if changes else "legacy_slot_match" if legacy else "match"
            ),
            differences=changes,
            observed_template=actual_template,
            observed_leaving={str(k): v for k, v in sorted(activated.leaving.items())},
        )
    except (ValueError, RuntimeError) as error:
        # The exact error text is useful to review but not a portable API contract.
        activation["error_type"] = type(error).__name__
    if activation["status"] in ("restoration_failed", "activation_failed", "mismatch"):
        issues.append(
            _issue(
                "activation_exception",
                "Automatic reactivation does not reproduce this stored definition.",
                status=activation["status"],
                differences=activation.get("differences", []),
            )
        )
    return {
        "definition_hash": definition_hash(mol),
        "issues": issues,
        "activation": activation,
    }


def audit_library(path=None):
    """Fresh full audit, always defaulting to the bundled rather than active SDF."""
    path = (
        Path(str(files("pyPept.data").joinpath("monomers.sdf")))
        if path is None
        else Path(path)
    )
    records = {}
    for index, molecule in enumerate(Chem.SDMolSupplier(str(path), removeHs=False)):
        if molecule is None:
            raise ValueError(f"Unreadable monomer at SDF record {index + 1}")
        symbol = _property(molecule, "symbol")
        if not symbol or symbol in records:
            raise ValueError(f"Missing or duplicated monomer symbol: {symbol!r}")
        records[symbol] = audit_monomer(molecule)
    records = dict(sorted(records.items()))
    return {
        "schema_version": 1,
        "audit_version": AUDIT_VERSION,
        "library_sha256": sha256(path.read_bytes()).hexdigest(),
        "generated_with_rdkit": rdBase.rdkitVersion,
        "scope": "Software compatibility facts, not scientific chemical validation.",
        "summary": {
            "records": len(records),
            "activation": dict(
                sorted(
                    Counter(
                        record["activation"]["status"] for record in records.values()
                    ).items()
                )
            ),
            "issues": dict(
                sorted(
                    Counter(
                        issue["code"]
                        for record in records.values()
                        for issue in record["issues"]
                    ).items()
                )
            ),
        },
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = audit_library()
    args.output.write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
