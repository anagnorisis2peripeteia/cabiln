"""CABILN conversion request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from pyPept.inputs import convert_input, format_source
from pyPept.smiles import convert_smiles

from .schemas import _ConvertReq, _SmilesToCabilnReq, _ToCabilnReq

router = APIRouter()


def _conversion_response(result):
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
    return response


@router.post("/convert_notation")
def convert_notation(req: _ConvertReq):
    try:
        notation = {"bracket": "bracket", "branch": "percent"}.get(req.target)
        if notation is None:
            return JSONResponse(
                {"error": "target must be 'bracket' or 'branch'"}, status_code=400
            )
        return {"result": format_source(req.cabiln, notation)}
    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


@router.post("/smiles_to_cabiln")
def smiles_to_cabiln_endpoint(req: _SmilesToCabilnReq):
    """Convert a peptide SMILES to CABILN notation using pyPept's monomer library."""
    try:
        return _conversion_response(convert_smiles(req.smiles, notation=req.notation))
    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


@router.post("/to_cabiln")
def to_cabiln_endpoint(req: _ToCabilnReq):
    """Convert SMILES, BILN, or HELM input to CABILN notation."""
    try:
        kind, result = convert_input(req.input, notation=req.notation)
        if kind == "SMILES":
            return {**_conversion_response(result), "from": kind}
        return {"cabiln": result, "from": kind, "details": []}
    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)
