"""CABILN rendering request handlers."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from .cache import _rc_get, _rc_put, library_version
from .drawing import _draw_mol, _mol_block
from .notation import parse_source
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

        from pyPept.molecule import Molecule

        messages = []
        seq, parsed_source = parse_source(req.cabiln, warning_sink=messages.append)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt="ROMol")
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

        crosslink_groups = _build_crosslink_groups(seq)
        bracket_groups = _build_bracket_groups(
            seq, chain_ids, parsed_source, crosslink_groups
        )

        result = {
            "svg": svg,
            "mol_block": block,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
            "residue_map": {str(k): v for k, v in res_map.items()},
            "residues": residues,
            "chains": chains,
            "bracket_groups": bracket_groups,
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
        from rdkit import Chem
        from rdkit.Chem.Descriptors import ExactMolWt

        romol = Chem.MolFromSmiles(req.smiles)
        if romol is None:
            return JSONResponse({"error": "Invalid SMILES"}, status_code=400)

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

    txt = req.input.strip()
    w, h = max(400, req.width), max(300, req.height)
    romol = None
    fmt = None

    romol = Chem.MolFromSmiles(txt)
    if romol is not None:
        fmt = "SMILES"

    if romol is None and "PEPTIDE" in txt.upper() and "$" in txt:
        try:

            from pyPept.converter import Converter
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence

            conv = Converter(helm=txt)
            biln = conv.get_biln()
            seq = Sequence(biln, fmt="biln")
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt="ROMol")
            fmt = "HELM"
        except Exception:
            pass

    if romol is None:
        try:
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence

            seq = Sequence(txt, fmt="biln")
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt="ROMol")
            fmt = "BILN"
        except Exception:
            pass

    if romol is None:
        try:
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence

            seq, _ = parse_source(txt)
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt="ROMol")
            fmt = "CABILN"
        except Exception:
            pass

    if romol is None:
        return JSONResponse(
            {"error": "Could not parse as SMILES, BILN, HELM, or CABILN"},
            status_code=400,
        )

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

        romol = Chem.MolFromMolBlock(req.mol_block)
        if romol is None:
            return JSONResponse({"error": "Invalid .mol data"}, status_code=400)

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
        from rdkit import Chem

        from pyPept.molecule import Molecule

        smiles_mol = Chem.MolFromSmiles(req.smiles)
        if smiles_mol is None:
            return JSONResponse({"error": "Invalid SMILES"}, status_code=400)

        seq, _ = parse_source(req.cabiln)
        mol = Molecule(seq)
        cabiln_mol = mol.get_molecule(fmt="ROMol")
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


def _build_bracket_groups(seq, chain_ids, cabiln="", crosslink_groups=None):
    """Identify bracket branch groups and their host monomers.

    Each separate attachment (``.[A.B.C]``, ``.[D]``, ``.cap(r,r)``) on a
    host residue becomes its own group.  Groups are split by connected
    components among the branch residues — branch residues that bond to
    each other stay together; those that only bond to the host get their
    own group.

    Returns list of ``{host, members}``.
    """
    if not chain_ids or len(chain_ids) < 2:
        return []
    main_set = set(chain_ids[0])
    branch_set = set()
    for ci in range(1, len(chain_ids)):
        branch_set.update(chain_ids[ci])
    if not branch_set:
        return []

    xlink_pairs = set()
    for g in crosslink_groups or []:
        if len(g["members"]) == 2:
            a, b = g["members"]
            xlink_pairs.add((a, b))
            xlink_pairs.add((b, a))

    host_of = {}
    branch_adj: dict[int, set] = {idx: set() for idx in branch_set}
    for bond in seq.s_bonds:
        m1, m2 = bond[0], bond[2]
        if (m1, m2) in xlink_pairs:
            continue
        if m1 in main_set and m2 in branch_set:
            host_of[m2] = m1
        elif m2 in main_set and m1 in branch_set:
            host_of[m1] = m2
        elif m1 in branch_set and m2 in branch_set:
            branch_adj[m1].add(m2)
            branch_adj[m2].add(m1)

    # Propagate host through inter-branch bonds
    changed = True
    while changed:
        changed = False
        for idx in branch_set:
            if idx in host_of:
                continue
            for peer in branch_adj.get(idx, set()):
                if peer in host_of:
                    host_of[idx] = host_of[peer]
                    changed = True
                    break

    # Connected components among branch residues (union-find)
    parent = {idx: idx for idx in branch_set}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for idx, peers in branch_adj.items():
        for p in peers:
            union(idx, p)

    from collections import defaultdict

    components = defaultdict(list)
    for idx in branch_set:
        if idx in host_of:
            components[(host_of[idx], find(idx))].append(idx)

    return [
        {"host": key[0], "members": sorted(members)}
        for key, members in components.items()
    ]


def _build_crosslink_groups(seq):
    """Extract crosslink !n pairs from the parsed BILN."""
    import re

    biln = seq.s_biln
    chains = biln.split(".")
    m_idx = 0
    tag_to_monomers = {}
    for chain in chains:
        residues = chain.split("-")
        for res in residues:
            for m in re.finditer(r"\((!\w+),\d+\)", res):
                tag = m.group(1)
                tag_to_monomers.setdefault(tag, []).append(m_idx)
            m_idx += 1
    return [
        {"tag": tag, "members": mems}
        for tag, mems in tag_to_monomers.items()
        if len(mems) == 2
    ]
