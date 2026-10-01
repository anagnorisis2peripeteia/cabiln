"""HTTP adapters for source-aware peptide editing and attachment inspection."""

from __future__ import annotations

from fastapi import APIRouter, Query
from fastapi.responses import JSONResponse

from pyPept.attachments import attachment_sites, reaction_for_types
from pyPept.editor import PeptideDocument
from pyPept.molecule import Molecule
from pyPept.monomer_store import _load_sdf
from pyPept.peptide import Peptide
from pyPept.sequence import Sequence

from .drawing import _draw_mol
from .cache import _rc_get, _rc_put, library_version
from .execution import error_response
from .projects import checked_context, project_context
from .schemas import (
    _InsertBackboneReq, _InsertBondReq, _ValidateBondReq,
    _ReplacementOptionsReq, _ReplaceMonomerReq,
)

router = APIRouter()


def _initial_monomer(symbol):
    sequence = Sequence(symbol)
    if len(sequence.s_monomers) != 1:
        raise ValueError("Choose a single monomer to insert")
    Molecule(sequence)
    return {"result": symbol}


@router.post("/insert_bond")
def insert_bond(req: _InsertBondReq):
    try:
        if not req.cabiln.strip():
            return _initial_monomer(req.new_abbr)
        document = PeptideDocument(req.cabiln)
        host = document.select(req.host_residue_idx)
        if req.target_residue_idx >= 0:
            result = document.connect(
                host, req.r_host, document.select(req.target_residue_idx), req.r_new
            )
        else:
            result = document.attach(host, req.r_host, req.new_abbr, req.r_new)
        return {"result": result}
    except Exception as exc:
        return error_response(exc)


@router.post("/insert_backbone")
def insert_backbone(req: _InsertBackboneReq):
    try:
        if not req.cabiln.strip():
            return _initial_monomer(req.new_abbr)
        document = PeptideDocument(req.cabiln)
        result = document.insert_backbone(document.select(req.after_idx), req.new_abbr)
        return {"result": result}
    except Exception as exc:
        return error_response(exc)


def _replacement_catalog():
    """Retain detected sites once per library revision, including concrete cap forms."""
    version = library_version()
    key = (version, "replacement_sites")
    cached = _rc_get(key)
    if cached is not None:
        return cached
    _, monomers = _load_sdf()
    catalog = []
    for abbr, molecule in monomers.items():
        props = molecule.GetPropsAsDict()
        sites = attachment_sites(molecule)
        catalog.append({
            "abbr": abbr, "name": props.get("m_name", ""),
            "type": props.get("m_type", ""), "subtype": props.get("m_subtype", ""),
            "chem_types": ",".join(f"{s['slot']}:{s['chem_type']}" for s in sites),
            "sites": [{"slot": s["slot"], "chem_type": s["chem_type"]} for s in sites],
        })
    if library_version() == version:
        _rc_put(key, catalog)
    return catalog


@router.post("/replacement_options")
def replacement_options(req: _ReplacementOptionsReq):
    try:
        context = project_context()
        document = PeptideDocument(req.cabiln)
        selection = document.select(req.residue_idx)
        requirements = document.replacement_requirements(selection)
        current = document.peptide.occurrence(req.residue_idx).symbol
        candidates = []
        for monomer in _replacement_catalog():
            if monomer["abbr"] == current:
                continue
            options = document.replacement_options(requirements, monomer["sites"])
            if options is not None:
                candidates.append({**monomer, **options})
        return {"source_echo": req.cabiln, "residue_idx": req.residue_idx,
                "requirements": requirements, "candidates": candidates,
                "context": checked_context(context)}
    except Exception as exc:
        return error_response(exc)


@router.post("/replace_monomer")
def replace_monomer(req: _ReplaceMonomerReq):
    try:
        context = checked_context(req.context)
        document = PeptideDocument(req.cabiln)
        result = document.replace_monomer(
            document.select(req.residue_idx), req.new_abbr, req.slot_map
        )
        return {"result": result, "source_echo": req.cabiln,
                "context": checked_context(context)}
    except Exception as exc:
        return error_response(exc)


@router.post("/validate_bond")
def validate_bond(req: _ValidateBondReq):
    """Check if a bond between two R-group chemistry types is valid."""
    try:
        entry = reaction_for_types(req.chem_type_a, req.chem_type_b)
        if entry:
            return {
                "valid": True,
                "reaction": entry.get("id", "unknown"),
                "description": entry.get("description", ""),
            }
        else:
            return {
                "valid": False,
                "reason": f"No reaction for {req.chem_type_a} + {req.chem_type_b}",
            }
    except Exception as exc:
        return error_response(exc)


@router.get("/reactions")
def list_reactions():
    """Return all valid (chem_type_a, chem_type_b) pairs from the reaction index."""
    try:
        from pyPept.interfaces.reaction_library import REACTION_INDEX

        pairs = [list(pair) for pair in REACTION_INDEX.keys()]
        return JSONResponse(pairs)
    except Exception as exc:
        return error_response(exc)


@router.get("/monomer_rgroups")
def monomer_rgroups(
    abbr: str = Query(max_length=100),
    residue_idx: int = Query(default=-1, ge=-1),
    cabiln: str = Query(default="", max_length=20000),
):
    """Inspect the selected instance, or a tile from the current monomer store."""
    try:
        from rdkit import Chem

        used_slots = set()
        leaving_groups = None
        if cabiln.strip():
            sequence = Sequence(cabiln)
            if not 0 <= residue_idx < len(sequence.s_monomers):
                raise ValueError(f"Residue {residue_idx} does not exist")
            monomer = sequence.s_monomers[residue_idx]
            target_mol = Chem.Mol(monomer["m_romol"])
            leaving_groups = monomer["m_Rgroups"]
            abbr = monomer["m_abbr"]
            used_slots = {
                endpoint.slot
                for endpoint in Peptide.occupied_sites_from_sequence(sequence)
                if endpoint.occurrence_id == residue_idx
            }
        else:
            if residue_idx >= 0:
                raise ValueError("A selected residue requires its sequence")
            _, monomers = _load_sdf()
            target_mol = monomers.get(abbr)
            if target_mol is None:
                return JSONResponse(
                    {"error": f"Monomer '{abbr}' not found"}, status_code=404
                )
            target_mol = Chem.Mol(target_mol)

        sites = attachment_sites(target_mol, leaving_groups)
        site_atoms = {
            atom.GetIsotope(): atom.GetIdx()
            for atom in target_mol.GetAtoms()
            if atom.GetAtomicNum() == 0
        }
        for site in sites:
            site["used"] = site["slot"] in used_slots
            site["atom_idx"] = site_atoms[site["slot"]]
        svg = _draw_mol(target_mol, 180, 140, used_slots or None)
        return {"svg": svg, "rgroups": sites, "abbr": abbr}
    except Exception as exc:
        return error_response(exc)
