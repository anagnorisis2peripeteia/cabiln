"""CABILN conversion request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from pyPept.smiles import convert_smiles

from .notation import format_source
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
    from rdkit import Chem as _C

    txt = req.input.strip()

    # 1. Try SMILES
    mol = _C.MolFromSmiles(txt)
    if mol is not None:
        try:
            return {
                **_conversion_response(convert_smiles(txt, notation=req.notation)),
                "from": "SMILES",
            }
        except Exception as exc:
            return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)

    # 2. Try HELM
    if "PEPTIDE" in txt.upper() and "$" in txt:
        try:

            from pyPept.converter import Converter
            from pyPept.sequence import biln_to_cabiln

            biln = Converter(helm=txt).get_biln()
            cabiln = biln_to_cabiln(biln)
            # Verify it assembles
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence

            Molecule(Sequence(cabiln)).get_molecule(fmt="ROMol")
            cabiln = format_source(cabiln, req.notation)
            return {"cabiln": cabiln, "from": "HELM", "details": []}
        except Exception as exc:
            return JSONResponse(
                {"error": f"HELM parse failed: {exc}".split("\n")[0]}, status_code=400
            )

    # 3. Try BILN
    try:
        from pyPept.molecule import Molecule
        from pyPept.sequence import Sequence, biln_to_cabiln

        cabiln = biln_to_cabiln(txt)
        Molecule(Sequence(cabiln)).get_molecule(fmt="ROMol")
        cabiln = format_source(cabiln, req.notation)
        return {"cabiln": cabiln, "from": "BILN", "details": []}
    except Exception:
        pass

    return JSONResponse(
        {"error": "Could not parse input as SMILES, BILN, or HELM"}, status_code=400
    )
