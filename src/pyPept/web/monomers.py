"""CABILN monomers request handlers."""

from __future__ import annotations

from hashlib import sha256

from fastapi import APIRouter, Query, Request, Response
from fastapi.responses import JSONResponse

from pyPept.attachments import attachment_sites
from pyPept.monomer_store import (
    _load_sdf,
    library_path,
    monomer_record,
    register_molecule,
)

from .cache import _rc_get, _rc_put, library_version
from .drawing import _draw_mol
from .execution import error_response
from .monomer_display import _restore_leaving_groups, _restore_reagent_form
from .schemas import _PreviewReq, _RegisterReq
from .security import require_registration

router = APIRouter()


@router.get("/monomers")
def list_monomers():
    try:
        return palette_response()
    except Exception as exc:
        return error_response(exc)


def palette_response():
    """Prepare JSON once per library revision; share the bounded cache with startup."""
    version = library_version()
    key = (version, "monomers")
    token = sha256(repr(version).encode("utf-8")).hexdigest()
    cached = _rc_get(key)
    if cached is not None:
        return Response(cached, media_type="application/json",
                        headers={"X-Library-Version": token})
    all_mols, mol_by_abbr = _load_sdf()
    from pyPept.sequence import get_monomer_info
    from pyPept.library_quality import quality_for_monomer

    aliases = get_monomer_info(str(library_path())).attrs.get("_degen_aliases", {})
    monomers = []
    degen_nterm = {}  # base -> entry with trailing _
    degen_cterm = {}  # base -> entry with leading _
    for mol in all_mols:
        if mol is None:
            continue
        p = mol.GetPropsAsDict()
        abbr = p.get("m_abbr", "") or p.get("symbol", "")
        if not abbr:
            continue
        rgroups = p.get("m_Rgroups", "")
        slots = [r.strip() for r in rgroups.split(",")]
        lg_parts = [
            f"R{i+1}:{slots[i]}"
            for i in range(len(slots))
            if slots[i] not in ("None", "", "none")
        ]
        sites = attachment_sites(mol)
        entry = {
            "abbr": abbr,
            "name": p.get("m_name", ""),
            "type": p.get("m_type", ""),
            "subtype": p.get("m_subtype", ""),
            "chem_types": ",".join(f"{s['slot']}:{s['chem_type']}" for s in sites),
            "declared_chem_types": p.get("m_chem_types", ""),
            "backbone_insertable": {1, 2}.issubset({s["slot"] for s in sites}),
            "leaving": ", ".join(lg_parts),
            "quality": quality_for_monomer(mol),
        }
        if abbr.endswith("_") and abbr[:-1] in aliases:
            degen_nterm[abbr[:-1]] = entry
            continue
        if abbr.startswith("_") and abbr[1:] in aliases:
            degen_cterm[abbr[1:]] = entry
            continue
        monomers.append(entry)
    all_bases = set(degen_nterm) | set(degen_cterm)
    for base in sorted(all_bases):
        nt = degen_nterm.get(base)
        ct = degen_cterm.get(base)
        primary = nt or ct
        merged = {
            "abbr": base,
            "name": primary["name"],
            "type": primary["type"],
            "subtype": primary["subtype"],
            "chem_types": ",".join(
                entry["chem_types"] for entry in (nt, ct) if entry
            ),
            "backbone_insertable": False,
            "leaving": primary["leaving"],
            "degenerate": True,
            "quality": _merged_quality(nt, ct),
        }
        if nt:
            merged["nterm_abbr"] = nt["abbr"]
            merged["nterm_leaving"] = nt["leaving"]
        if ct:
            merged["cterm_abbr"] = ct["abbr"]
            merged["cterm_leaving"] = ct["leaving"]
        monomers.append(merged)
    # A definition can change without changing its tile labels. Expose its
    # revision so browser previews can invalidate independently of the body.
    # Do not certify or cache a response assembled across a library update.
    response = JSONResponse(monomers)
    if library_version() == version:
        _rc_put(key, response.body)
        response.headers["X-Library-Version"] = token
    return response


def _merged_quality(*entries):
    quality = [entry for entry in entries if entry]
    statuses = {entry["quality"]["status"] for entry in quality}
    status = (
        "unreviewed"
        if "unreviewed" in statuses
        else (
            "review_required" if "review_required" in statuses else "no_known_exception"
        )
    )
    return {
        "status": status,
        "issues": [
            {**issue, "message": f"{entry['abbr']}: {issue['message']}"}
            for entry in quality
            for issue in entry["quality"]["issues"]
        ],
    }


@router.get("/monomer_svg")
def monomer_svg(
    abbr: str = Query(max_length=100),
    width: int = Query(default=220, ge=64, le=4096),
    height: int = Query(default=180, ge=64, le=4096),
):
    _mck = (library_version(), "monomer", abbr, width, height)
    _mhit = _rc_get(_mck)
    if _mhit is not None:
        return _mhit
    try:
        from rdkit import Chem

        _all_mols, mol_by_abbr = _load_sdf()

        if abbr in mol_by_abbr:
            mol = Chem.Mol(mol_by_abbr[abbr])
            svg = _draw_mol(mol, width, height)
            restored = _restore_leaving_groups(mol)
            svg_restored = _draw_mol(restored, width, height)
            result = {"svg": svg, "svg_restored": svg_restored}

            reagent_mol, reagent_meta = _restore_reagent_form(mol, abbr)
            if reagent_mol is not None:
                result["svg_reagent"] = _draw_mol(reagent_mol, width, height)
                result["reagent"] = reagent_meta
            _rc_put(_mck, result)
            return result

        # Check if this is a degenerate base name (e.g. "Bn" -> Bn_/_Bn)
        # Detect pairs: look for abbr_ and _abbr variants
        nterm_key = abbr + "_"
        cterm_key = "_" + abbr
        variants_found = []
        if nterm_key in mol_by_abbr:
            variants_found.append(nterm_key)
        if cterm_key in mol_by_abbr:
            variants_found.append(cterm_key)

        if len(variants_found) >= 2:
            panels = []
            for vkey in variants_found:
                vmol = mol_by_abbr[vkey]
                restored_v = _restore_leaving_groups(vmol)
                label = f"N-term ({vkey})" if vkey.endswith("_") else f"C-term ({vkey})"
                reagent_mol, reagent_meta = _restore_reagent_form(vmol, vkey)
                panel = {
                    "label": label,
                    "svg": _draw_mol(restored_v, width, height),
                }
                if reagent_mol is not None:
                    panel["svg_reagent"] = _draw_mol(reagent_mol, width, height)
                    panel["reagent"] = reagent_meta
                panels.append(panel)
            rgroup_svg = _draw_mol(
                Chem.Mol(mol_by_abbr[variants_found[0]]), width, height
            )
            return {
                "degenerate": True,
                "variants": panels,
                "svg": rgroup_svg,
            }

        return JSONResponse({"error": f"Monomer '{abbr}' not found"}, status_code=404)
    except Exception as exc:
        return error_response(exc)


@router.post("/preview_monomer")
def preview_monomer(req: _PreviewReq):
    try:

        from pyPept.interfaces.monomer_pipeline import pre_activate

        result = pre_activate(req.smiles)

        from rdkit import Chem

        mol = Chem.MolFromSmiles(result.chuckles)
        if mol is None:
            return JSONResponse(
                {"error": "Generated CHUCKLES is invalid"}, status_code=400
            )

        svg = _draw_mol(mol, req.width, req.height)
        return {
            "chuckles": result.chuckles,
            "chem_types": {str(k): v for k, v in result.chem_types.items()},
            "leaving": {str(k): v for k, v in result.leaving.items()},
            "svg": svg,
        }

    except Exception as exc:
        if isinstance(exc, ValueError) and str(exc).startswith(
            "Cannot identify backbone (N + COOH) or a single cap terminus."
        ):
            return JSONResponse(
                {"error": (
                    "No supported attachment sites were detected. Preview a "
                    "monomer with a supported backbone, cap or sidechain site."
                )},
                status_code=400,
            )
        return error_response(exc)


@router.post("/register_monomer")
def register_monomer(req: _RegisterReq, request: Request):
    require_registration(request)
    try:
        from rdkit.Chem import rdDepictor

        # The preview payload describes every detected site. Bulk import also
        # accepts standalone records and sparse declarations from older SDFs.
        if not req.leaving:
            raise ValueError(
                "Each attachment needs a unique numbered dummy with one neighbour"
            )
        if set(req.chem_types) != set(req.leaving):
            raise ValueError(
                "Attachment slots, chemistry types, and leaving groups must match"
            )
        rdDepictor.SetPreferCoordGen(True)
        mol = monomer_record(
            req.chuckles,
            req.abbr,
            req.leaving,
            req.chem_types,
            name=req.name,
            m_type=req.type,
            m_subtype=req.subtype,
            minimum_slots=6,
        )

        return {"ok": True, "total": register_molecule(mol)}

    except Exception as exc:
        return error_response(exc)
