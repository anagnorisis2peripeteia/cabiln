"""CABILN conversion request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from pyPept.inputs import convert_input, format_source
from pyPept.smiles import convert_smiles

from .execution import error_response
from .projects import checked_context, project_context
from .schemas import _ConvertReq, _SmilesToCabilnReq, _ToCabilnReq

router = APIRouter()


def _conversion_response(result, context=None):
    response = {
        "cabiln": result.cabiln,
        "details": [{"abbr": a, "score": s, "total": t} for a, s, t in result.details],
        "inferred_stereo": result.inferred_stereo,
        "synthetic_components": list(result.synthetic_components),
        "warnings": list(result.warnings),
        "recognition_status": result.recognition_status,
        "search_complete": result.search_complete,
        "assignments": [
            {
                "symbol": assignment.symbol,
                "source_atoms": list(assignment.source_atoms),
                "attachments": list(assignment.attachments),
                "recognized": assignment.recognized,
                "residue_index": assignment.residue_index,
            }
            for assignment in result.assignments
        ],
    }
    if result.warnings:
        response["warning"] = " ".join(result.warnings)
    if context is not None:
        response["context"] = checked_context(context)
    return response


@router.post("/convert_notation")
def convert_notation(req: _ConvertReq):
    try:
        notation = {"bracket": "bracket", "branch": "percent"}.get(req.target)
        if notation is None:
            return JSONResponse(
                {"error": "target must be 'bracket' or 'branch'"}, status_code=400
            )
        context = project_context()
        result = format_source(req.cabiln, notation, canonical=req.canonical)
        response = {"result": result, "context": checked_context(context)}
        if req.canonical:
            from pyPept.canonical import canonical_convention

            response["canonical"] = {
                **canonical_convention(),
                "binding": context["library_binding"],
            }
        return response
    except Exception as exc:
        return error_response(exc)


@router.post("/smiles_to_cabiln")
def smiles_to_cabiln_endpoint(req: _SmilesToCabilnReq):
    """Convert a peptide SMILES to CABILN notation using pyPept's monomer library."""
    try:
        context = project_context()
        return _conversion_response(
            convert_smiles(req.smiles, notation=req.notation), context
        )
    except Exception as exc:
        return error_response(exc)


@router.post("/to_cabiln")
def to_cabiln_endpoint(req: _ToCabilnReq):
    """Convert the selected or detected input format to CABILN notation."""
    try:
        context = project_context()
        kind, result = convert_input(
            req.input, notation=req.notation, input_format=req.input_format
        )
        if kind == "SMILES":
            return {**_conversion_response(result, context), "from": kind}
        return {
            "cabiln": result,
            "from": kind,
            "details": [],
            "context": checked_context(context),
        }
    except Exception as exc:
        return error_response(exc)
