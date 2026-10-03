"""Portable document validation bound to the definitions actually resolved."""

from __future__ import annotations

import json
from copy import deepcopy
from hashlib import sha256

from fastapi import APIRouter
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from pyPept.canonical import canonical_convention
from pyPept.monomer_store import library_binding

router = APIRouter()
FORMATS = {"cabiln", "smiles", "biln", "helm"}
MAX_PROJECT_BYTES = 2 * 1024 * 1024


class ProjectRequest(BaseModel):
    project: dict


class ProjectError(ValueError):
    def __init__(self, message, code="project_invalid"):
        super().__init__(message)
        self.code = code


def project_context():
    return {
        "project_version": 1,
        "library_binding": library_binding(),
        "canonical": canonical_convention(),
    }


def checked_context(context):
    if project_context() != context:
        raise ValueError("The monomer library changed; retry the operation")
    return context


def _text(value, field, maximum=20000):
    if not isinstance(value, str) or len(value) > maximum:
        raise ProjectError(f"{field} must be text of at most {maximum:,} characters")
    return value


def _document(value, field):
    if not isinstance(value, dict):
        raise ProjectError(f"{field} must be a document")
    _text(value.get("text"), f"{field}.text")
    if value.get("notation") not in FORMATS:
        raise ProjectError(f"{field}.notation is not supported")
    _text(value.get("warning", ""), f"{field}.warning")
    for key in ("quality", "canonical", "context"):
        if value.get(key) is not None and not isinstance(value[key], dict):
            raise ProjectError(f"{field}.{key} must be an object or null")
    quality = value.get("quality")
    if quality is not None:
        if quality.get("recognition_status") not in {
            "complete",
            "partial",
            "unresolved",
        }:
            raise ProjectError(f"{field}.quality has an invalid recognition status")
        for key in ("search_complete", "inferred_stereo"):
            if type(quality.get(key)) is not bool:
                raise ProjectError(f"{field}.quality.{key} must be a boolean")
        assignments = quality.get("assignments", [])
        if not isinstance(assignments, list) or len(assignments) > 20000:
            raise ProjectError(f"{field}.quality has invalid assignments")
    return value


def _shape(project):
    pending = [(project, 0)]
    while pending:
        value, depth = pending.pop()
        if depth > 32:
            raise ProjectError("Project metadata is nested too deeply")
        if isinstance(value, dict):
            pending.extend((item, depth + 1) for item in value.values())
        elif isinstance(value, list):
            pending.extend((item, depth + 1) for item in value)
    try:
        serialized = json.dumps(project, allow_nan=False).encode("utf-8")
    except (ValueError, TypeError, RecursionError) as exc:
        raise ProjectError("Project must contain finite JSON data") from exc
    if len(serialized) > MAX_PROJECT_BYTES:
        raise ProjectError("Project exceeds the 2 MiB limit")
    if (
        project.get("format") != "cabiln-project"
        or type(project.get("version")) is not int
        or project["version"] != 1
    ):
        raise ProjectError("Unsupported project format or version")
    project = deepcopy(project)
    documents = {"document": _document(project.get("document"), "document")}
    drafts = project.get("drafts", {})
    if not isinstance(drafts, dict) or set(drafts) - FORMATS:
        raise ProjectError("Project has unsupported notation drafts")
    for notation, draft in drafts.items():
        if isinstance(draft, dict):
            draft.setdefault("notation", notation)
        documents[f"draft:{notation}"] = _document(draft, f"drafts.{notation}")
        if draft["notation"] != notation:
            raise ProjectError("Draft notation does not match its key")
    reference = project.get("reference", {"text": "", "original": None})
    if not isinstance(reference, dict):
        raise ProjectError("Project reference must be an object")
    if reference.get("context") is not None and not isinstance(
        reference["context"], dict
    ):
        raise ProjectError("Project reference context must be an object or null")
    project["reference"] = reference
    _text(reference.get("text", ""), "reference.text")
    original = reference.get("original")
    if original is not None:
        if not isinstance(original, dict) or original.get("kind") not in {
            "text",
            "mol",
        }:
            raise ProjectError("Invalid original reference")
        _text(original.get("content"), "reference.original.content", 1000000)
        _text(original.get("name", ""), "reference.original.name", 256)
    _text(project.get("saved_at", ""), "saved_at", 100)
    return project, documents, reference


def _digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return sha256(data.encode("utf-8")).hexdigest()


def _resolution(document):
    from rdkit import Chem

    from pyPept.inputs import read_input
    from pyPept.peptide import Peptide
    from pyPept.structure import require_supported_stereo

    source = (document["notation"], document["text"])
    result = {"source": _digest(source), "resolved": None}
    if not document["text"].strip():
        result["resolved"] = _digest("empty")
        return result
    try:
        parsed = read_input(
            document["text"],
            input_format=document["notation"],
            warning_sink=lambda _: None,
        )
        if parsed.sequence is None:
            require_supported_stereo(parsed.molecule)
            identity = Chem.MolToCXSmiles(parsed.molecule)
        else:
            peptide = Peptide.from_sequence(parsed.sequence)
            identity = {
                "definitions": [occ.definition.key for occ in peptide.occurrences],
                "connections": [
                    [(end.occurrence_id, end.slot) for end in edge.endpoints]
                    for edge in peptide.connections
                ],
            }
        result["resolved"] = _digest(identity)
    except ValueError:
        # Invalid drafts are work worth saving, but cannot be rebound safely.
        pass
    return result


def _reference_resolution(reference):
    from rdkit import Chem

    from pyPept.inputs import detect_input, read_input
    from pyPept.structure import require_supported_stereo

    original = reference.get("original")
    text = original["content"] if original else reference.get("text", "")
    kind = original["kind"] if original else "text"
    result = {
        "source": _digest((kind, text, reference.get("text", ""))),
        "resolved": None,
    }
    if not text.strip():
        result["resolved"] = _digest("empty")
        return result
    try:
        parsed = (
            read_input(text, input_format="mol")
            if kind == "mol"
            else detect_input(text)
        )
        molecule = parsed.assemble(depiction="none")
        require_supported_stereo(molecule)
        result["resolved"] = _digest(Chem.MolToCXSmiles(molecule))
    except ValueError:
        pass
    return result


def _compatible(saved, current, signatures, name=None, *, preparing=False):
    if not isinstance(saved, dict) or saved.get("project_version") != 1:
        raise ProjectError("Project has no library context; reopen the original source")
    if saved.get("canonical") != current["canonical"]:
        raise ProjectError(
            "Project uses a different canonical/RDKit convention. Open its source "
            "as a new document and verify it before saving in this environment.",
            "project_binding_mismatch",
        )
    binding = saved.get("library_binding")
    saved_signatures = saved.get("resolutions", {})
    if not isinstance(saved_signatures, dict):
        raise ProjectError("Project contains invalid resolution metadata")
    selected = {name: signatures[name]} if name else signatures
    if binding == current["library_binding"]:
        if not preparing and any(
            saved_signatures.get(key) != signature
            for key, signature in selected.items()
        ):
            raise ProjectError(
                "Project sources differ from their saved chemistry context. Open "
                "the edited source as a new document and verify it before saving.",
                "project_binding_mismatch",
            )
        return
    if not isinstance(binding, dict) or any(
        binding.get(key) != current["library_binding"][key]
        for key in ("reactions", "caps", "site_chemistry")
    ):
        raise ProjectError(
            "Project chemistry rules differ from this installation. Use its original "
            "library/release. If attachment numbers changed, migrate its notation "
            "with pyPept.monomer_migration before opening it here.",
            "project_binding_mismatch",
        )
    for key, signature in selected.items():
        if signature["resolved"] is None or saved_signatures.get(key) != signature:
            raise ProjectError(
                "Project library definitions differ or a draft cannot be resolved. "
                "Use the original library, or migrate its notation with "
                "pyPept.monomer_migration before saving a new project. "
                "No current document was replaced.",
                "project_binding_mismatch",
            )


def validate_project(project, *, preparing=False):
    project, documents, reference = _shape(project)
    current = project_context()
    signatures = {name: _resolution(doc) for name, doc in documents.items()}
    signatures["reference"] = _reference_resolution(reference)
    saved_context = project.get("context")
    if saved_context is not None or not preparing:
        _compatible(saved_context, current, signatures, preparing=preparing)
    for name, document in documents.items():
        context = document.get("context")
        if context is not None:
            resolutions = context.get("resolutions", {})
            if not isinstance(resolutions, dict):
                raise ProjectError("Document contains invalid resolution metadata")
            previous = resolutions.get(name)
            if previous is not None and not isinstance(previous, dict):
                raise ProjectError("Document contains an invalid resolution")
            if (
                preparing
                and previous is not None
                and previous.get("source") != signatures[name]["source"]
                and (document.get("quality") or document.get("canonical"))
            ):
                raise ProjectError(
                    "Document text changed since its chemistry metadata was saved. "
                    "Render and verify the edited source before saving.",
                    "project_binding_mismatch",
                )
            # Render contexts have no resolutions; the saved project context does.
            evidence = dict(context)
            if (
                isinstance(saved_context, dict)
                and (
                    saved_context.get("library_binding")
                    == context.get("library_binding")
                )
                and saved_context.get("resolutions")
            ):
                evidence["resolutions"] = saved_context["resolutions"]
            _compatible(evidence, current, signatures, name, preparing=preparing)
        elif preparing and (document.get("quality") or document.get("canonical")):
            raise ProjectError(
                "Render this document before saving its chemistry metadata.",
                "project_binding_mismatch",
            )
    reference_context = reference.get("context")
    if reference_context is not None:
        evidence = dict(reference_context)
        if (
            isinstance(saved_context, dict)
            and saved_context.get("library_binding")
            == reference_context.get("library_binding")
            and saved_context.get("resolutions")
        ):
            evidence["resolutions"] = saved_context["resolutions"]
        _compatible(evidence, current, signatures, "reference", preparing=preparing)
    if project_context() != current:
        raise ProjectError("Library changed while checking the project; retry saving")
    current["resolutions"] = signatures
    project["context"] = current
    for document in documents.values():
        document["context"] = current
    reference["context"] = current
    if len(json.dumps(project, allow_nan=False).encode("utf-8")) > MAX_PROJECT_BYTES:
        raise ProjectError("Prepared project exceeds the 2 MiB limit")
    return project


def _project_error(exc):
    if isinstance(exc, ProjectError):
        return JSONResponse({"error": str(exc), "code": exc.code}, status_code=409)
    from .execution import error_response

    return error_response(exc)


@router.get("/project_context")
def get_project_context():
    try:
        return project_context()
    except Exception as exc:
        return _project_error(exc)


@router.post("/prepare_project")
def prepare_project(req: ProjectRequest):
    try:
        return {"project": validate_project(req.project, preparing=True)}
    except Exception as exc:
        return _project_error(exc)


@router.post("/validate_project")
def check_project(req: ProjectRequest):
    try:
        project = validate_project(req.project)
        return {"valid": True, "context": project["context"]}
    except Exception as exc:
        return _project_error(exc)
