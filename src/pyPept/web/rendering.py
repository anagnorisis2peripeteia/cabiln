"""CABILN rendering request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from pyPept.inputs import detect_input, read_input
from pyPept.structure import require_supported_stereo

from .cache import _rc_get, _rc_put, library_version
from .drawing import _draw_mol, _mol_block
from .execution import error_response
from .projects import checked_context, project_context
from .schemas import _CabilnReq, _MolBlockReq, _ReferenceReq, _SmilesReq, _VerifyReq

router = APIRouter()


@router.post("/render")
def render(req: _CabilnReq):
    w, h = max(400, req.width), max(300, req.height)
    try:
        from rdkit.Chem.Descriptors import ExactMolWt

        from pyPept.canonical import canonical_convention
        from pyPept.monomer_store import library_binding
        from pyPept.peptide import Peptide

        binding = library_binding()
        _ck = (library_version(), req.cabiln, w, h, req.seed)
        _hit = _rc_get(_ck)
        if (
            _hit is not None
            and _hit.get("context", {}).get("library_binding") == binding
        ):
            if library_binding() != binding:
                raise ValueError("The monomer library changed; retry the render")
            return _hit
        messages = []
        parsed = read_input(req.cabiln, warning_sink=messages.append, track_source=True)
        seq, parsed_source = parsed.sequence, parsed.source
        romol = parsed.assemble(depiction=None)
        mol = parsed.assembly
        if romol is None:
            return JSONResponse(
                {"error": "Assembly produced no molecule"}, status_code=400
            )

        svg = _draw_mol(romol, w, h, seed=req.seed)
        block = _mol_block(romol)

        res_map = mol.get_residue_atom_map()
        peptide = Peptide.from_sequence(seq)

        result = {
            "svg": svg,
            "mol_block": block,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
            "residue_map": {str(k): v for k, v in res_map.items()},
            **_render_presentation(seq, peptide, messages),
            "context": {
                "project_version": 1,
                "library_binding": binding,
                "canonical": canonical_convention(),
            },
        }
        if parsed_source != req.cabiln:
            result["normalized_cabiln"] = parsed_source
            result["normalization_note"] = "Converted legacy positional notation to bracket form."
        if library_binding() != binding:
            raise ValueError("The monomer library changed; retry the render")
        _rc_put(_ck, result)
        return result

    except Exception as exc:
        return error_response(exc)


def _structure_view(romol, req):
    from rdkit.Chem.Descriptors import ExactMolWt

    require_supported_stereo(romol)
    return {
        "svg": _draw_mol(romol, max(400, req.width), max(300, req.height)),
        "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
    }


@router.post("/render_smiles")
def render_smiles(req: _SmilesReq):
    try:
        context = project_context()
        romol = read_input(req.smiles, input_format="smiles").molecule
        return {
            **_structure_view(romol, req),
            "context": checked_context(context),
        }

    except Exception as exc:
        return error_response(exc)


@router.post("/render_reference")
def render_reference(req: _ReferenceReq):
    """Render the requested input format, using auto-detection by default."""
    try:
        from rdkit import Chem

        context = project_context()
        parsed = (
            detect_input(req.input, depiction=None)
            if req.input_format == "auto"
            else read_input(req.input, input_format=req.input_format)
        )
        romol, fmt = parsed.assemble(depiction=None), parsed.format
        return {
            **_structure_view(romol, req),
            "smiles": Chem.MolToSmiles(romol),
            "format": fmt,
            "context": checked_context(context),
        }
    except Exception as exc:
        return error_response(exc)


@router.post("/render_mol")
def render_mol(req: _MolBlockReq):
    try:
        from rdkit import Chem

        context = project_context()
        romol = read_input(req.mol_block, input_format="mol").molecule
        return {
            **_structure_view(romol, req),
            "smiles": Chem.MolToSmiles(romol),
            "context": checked_context(context),
        }

    except Exception as exc:
        return error_response(exc)


@router.post("/verify")
def verify(req: _VerifyReq):
    try:
        smiles_mol = read_input(req.smiles, input_format="smiles").molecule
        cabiln_mol = read_input(req.cabiln).assemble(depiction=None)
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
        return error_response(exc)


def _render_presentation(sequence, peptide, warnings):
    """Describe target parser IDs consistently for renders and conversions."""
    layout, crosslink_groups = _renderer_layout(peptide)
    return {
        "residues": [
            {
                "idx": i,
                "abbr": monomer.get("m_abbr", f"?{i}"),
                "kind": (
                    "synthetic" if occurrence.definition.synthetic_role else "library"
                ),
            }
            for i, (monomer, occurrence) in enumerate(
                zip(sequence.s_monomers, peptide.occurrences)
            )
        ],
        "chains": [
            {"idx": i, "residues": ids}
            for i, ids in enumerate(sequence.s_chains.get("s_monomerIDs", []))
        ],
        "layout": layout,
        "bracket_groups": layout["groups"],
        "crosslink_groups": crosslink_groups,
        "warnings": list(warnings),
        "cabiln_echo": str(sequence.s_inputbiln),
    }


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
