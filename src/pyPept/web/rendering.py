"""CABILN rendering request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from pyPept.inputs import detect_input, read_input

from .cache import _rc_get, _rc_put, library_version
from .drawing import _draw_mol, _mol_block
from .schemas import _CabilnReq, _MolBlockReq, _ReferenceReq, _SmilesReq, _VerifyReq

router = APIRouter()


@router.post("/render")
def render(req: _CabilnReq):
    w, h = max(400, req.width), max(300, req.height)
    _ck = (library_version(), req.cabiln, w, h, req.seed)
    _hit = _rc_get(_ck)
    if _hit is not None:
        return _hit
    try:
        from rdkit.Chem.Descriptors import ExactMolWt

        from pyPept.peptide import Peptide

        messages = []
        parsed = read_input(req.cabiln, warning_sink=messages.append, track_source=True)
        seq, parsed_source = parsed.sequence, parsed.source
        romol = parsed.assemble()
        mol = parsed.assembly
        if romol is None:
            return JSONResponse(
                {"error": "Assembly produced no molecule"}, status_code=400
            )

        svg = _draw_mol(romol, w, h, seed=req.seed)
        block = _mol_block(romol)

        res_map = mol.get_residue_atom_map()
        residues = [
            {"idx": i, "abbr": m.get("m_abbr", f"?{i}")}
            for i, m in enumerate(seq.s_monomers)
        ]

        chain_ids = seq.s_chains.get("s_monomerIDs", [])
        chains = [{"idx": ci, "residues": ids} for ci, ids in enumerate(chain_ids)]

        layout, crosslink_groups = _renderer_layout(Peptide.from_sequence(seq))

        result = {
            "svg": svg,
            "mol_block": block,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
            "residue_map": {str(k): v for k, v in res_map.items()},
            "residues": residues,
            "chains": chains,
            "layout": layout,
            "bracket_groups": layout["groups"],
            "crosslink_groups": crosslink_groups,
            "warnings": messages,
            "cabiln_echo": parsed_source,
        }
        if parsed_source != req.cabiln:
            result["normalized_cabiln"] = parsed_source
        _rc_put(_ck, result)
        return result

    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


@router.post("/render_smiles")
def render_smiles(req: _SmilesReq):
    try:
        from rdkit.Chem.Descriptors import ExactMolWt

        romol = read_input(req.smiles, input_format="smiles").molecule

        w, h = max(400, req.width), max(300, req.height)
        svg = _draw_mol(romol, w, h)
        return {
            "svg": svg,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
        }

    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


@router.post("/render_reference")
def render_reference(req: _ReferenceReq):
    """Auto-detect input format (SMILES, old BILN, HELM, CABILN) and render."""
    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    w, h = max(400, req.width), max(300, req.height)
    try:
        parsed = detect_input(req.input)
    except ValueError as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)
    romol, fmt = parsed.molecule, parsed.format

    svg = _draw_mol(romol, w, h)
    smiles = Chem.MolToSmiles(romol)
    return {
        "svg": svg,
        "smiles": smiles,
        "format": fmt,
        "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
    }


@router.post("/render_mol")
def render_mol(req: _MolBlockReq):
    try:
        from rdkit import Chem
        from rdkit.Chem.Descriptors import ExactMolWt

        romol = read_input(req.mol_block, input_format="mol").molecule

        w, h = max(400, req.width), max(300, req.height)
        svg = _draw_mol(romol, w, h)
        smiles = Chem.MolToSmiles(romol)
        return {
            "svg": svg,
            "smiles": smiles,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
        }

    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


@router.post("/verify")
def verify(req: _VerifyReq):
    try:
        smiles_mol = read_input(req.smiles, input_format="smiles").molecule
        cabiln_mol = read_input(req.cabiln).assemble()
        if cabiln_mol is None:
            return JSONResponse(
                {"error": "CABILN assembly produced no molecule"}, status_code=400
            )

        from pyPept.structure import compare_structures

        comparison = compare_structures(smiles_mol, cabiln_mol)

        resp = {
            "match": comparison.exact,
            "stereo_compatible": comparison.compatible,
            "smiles_canonical": comparison.source_smiles,
            "cabiln_canonical": comparison.candidate_smiles,
        }
        if comparison.compatible and not comparison.exact:
            resp["warning"] = (
                "The reference leaves some stereochemistry unspecified. "
                "Connectivity and specified stereochemistry agree; "
                "this is not an exact structure match."
            )
        return resp

    except Exception as exc:
        return JSONResponse({"error": str(exc).split("\n")[0]}, status_code=400)


def _renderer_layout(peptide):
    """Project the recorded source layout onto selectable residue IDs.

    Assembly chains describe chemical connectivity. Source segments and groups
    describe where occurrences appear, including nested and protected branches.
    """
    layout = peptide.layout
    crosslinks = [
        {
            "tag": edge.label,
            "members": list(dict.fromkeys(e.occurrence_id for e in edge.endpoints)),
        }
        for edge in peptide.connections
        if edge.label is not None
    ]
    links_by_tag = {item["tag"]: item["members"] for item in crosslinks}
    markers = [
        {
            "residue": marker.endpoint.occurrence_id,
            "tag": marker.label,
            "members": links_by_tag[marker.label],
            "group": marker.group,
            "before": marker.span.end
            <= peptide.occurrence(marker.endpoint.occurrence_id).source.token.start,
        }
        for marker in sorted(layout.markers, key=lambda marker: marker.span.start)
    ]

    def source_order(identities):
        return sorted(
            identities,
            key=lambda identity: peptide.occurrence(identity).source.token.start,
        )

    groups = []
    for group in layout.groups:
        nested_members = {
            member
            for child in layout.groups
            if child.parent == group.id
            for member in child.members
        }
        groups.append(
            {
                "id": group.id,
                "host": group.host,
                "members": list(group.members),
                "roots": source_order(set(group.members) - nested_members),
                "parent": group.parent,
                "kind": group.kind,
                "opening": group.opening,
                "closing": group.closing,
                "protected": group.protected,
            }
        )
    return {
        "segments": [
            {
                "id": segment.id,
                "roots": source_order(segment.roots),
                "members": list(segment.members),
            }
            for segment in layout.segments
        ],
        "groups": groups,
        "markers": markers,
    }, crosslinks
