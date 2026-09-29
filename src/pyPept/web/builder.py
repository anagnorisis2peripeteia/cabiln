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
from .execution import error_response
from .schemas import _InsertBackboneReq, _InsertBondReq, _ValidateBondReq

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
