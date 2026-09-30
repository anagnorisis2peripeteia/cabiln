"""Molecule depiction and comparison helpers."""

from __future__ import annotations


def _tag_atom_joins(svg, drawer, atom_count):
    """Give RDKit's untagged bond-join paths the identity of their exact vertex."""
    import re

    # RDKit smoothBondJoins draws a three-point path around an unlabelled atom.
    # Its middle point is that atom, rounded to one decimal place in the SVG.
    # Coincident positions are ambiguous, so do not assign those paths an owner.
    positions = {}
    for index in range(atom_count):
        point = drawer.GetDrawCoords(index)
        position = (round(point.x, 1), round(point.y, 1))
        positions[position] = index if position not in positions else None

    number = r"-?\d+\.\d+"
    join = (
        rf"<path d='M {number},{number} L ({number}),({number})"
        rf" L {number},{number}' style='fill:none;[^']*stroke-miterlimit:10;[^']*' />"
    )

    def tag(match):
        index = positions.get((float(match[1]), float(match[2])))
        if index is None:
            return match[0]
        return match[0].replace("<path ", f"<path class='atom-{index}' ", 1)

    return re.sub(join, tag, svg)


def _tag_stereo_annotations(svg, molecule):
    """Attach stereo glyphs to their atoms using RDKit's annotation draw order."""
    import re

    from rdkit import Chem

    marker = "class='CIP_Code'"
    if marker not in svg:
        return svg
    relative_groups = [
        group for group in molecule.GetStereoGroups()
        if group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE
    ]
    hidden_atoms = {atom.GetIdx() for group in relative_groups for atom in group.GetAtoms()}
    hidden_bonds = {bond.GetIdx() for group in relative_groups for bond in group.GetBonds()}
    owners = []
    # RDKit extractCIPCodes emits atom annotations, then bond annotations, in
    # index order. Both SVG text backends emit one element per character.
    for atom in molecule.GetAtoms():
        if atom.HasProp("_CIPCode") and atom.GetIdx() not in hidden_atoms:
            owners.extend([f"atom-{atom.GetIdx()}"] * (len(atom.GetProp("_CIPCode")) + 2))
    for bond in molecule.GetBonds():
        if bond.GetIdx() in hidden_bonds:
            continue
        code = bond.GetProp("_CIPCode") if bond.HasProp("_CIPCode") else {
            Chem.BondStereo.STEREOE: "E", Chem.BondStereo.STEREOZ: "Z",
        }.get(bond.GetStereo(), "")
        if code:
            owner = f"bond-{bond.GetIdx()} atom-{bond.GetBeginAtomIdx()} atom-{bond.GetEndAtomIdx()}"
            owners.extend([owner] * (len(code) + 2))
    # An unfamiliar RDKit output must not silently shift annotation ownership.
    # Generic atomNote/bondNote elements have class 'note' and remain untouched.
    if svg.count(marker) != len(owners):
        return svg
    owners = iter(owners)
    return re.sub(marker, lambda _: f"class='CIP_Code {next(owners)}'", svg)


def _indigo_layout(romol):
    """Lay out romol using Indigo's algorithm; copy coords back preserving atom indices."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.rdchem import Conformer

    try:
        from indigo import Indigo as _Indigo
    except ImportError:
        return None
    # mol block needs a conformer — compute a quick one just for the connectivity export
    tmp = Chem.RWMol(romol)
    if tmp.GetNumConformers() == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(tmp)
    mol_block = Chem.MolToMolBlock(tmp)
    try:
        indigo = _Indigo()
        indigo.setOption("ignore-stereochemistry-errors", True)
        im = indigo.loadMolecule(mol_block)
        im.layout()
        result_block = im.molfile()
        result = Chem.MolFromMolBlock(result_block, removeHs=False, sanitize=False)
        if result is None or result.GetNumConformers() == 0:
            return None
        src = result.GetConformer()
        n = romol.GetNumAtoms()
        if romol.GetNumConformers() == 0:
            conf = Conformer(n)
            for i in range(n):
                p = src.GetAtomPosition(i)
                conf.SetAtomPosition(i, (p.x, p.y, 0.0))
            romol.AddConformer(conf, assignId=True)
        else:
            conf = romol.GetConformer()
            for i in range(n):
                p = src.GetAtomPosition(i)
                conf.SetAtomPosition(i, (p.x, p.y, 0.0))
        return romol
    except Exception:
        return None


def _overlap_score(romol) -> int:
    """Count non-bonded atom pairs closer than 0.8 Å — lower is better."""
    conf = romol.GetConformer()
    n = romol.GetNumAtoms()
    pos = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(n)]
    bonded: set = set()
    for b in romol.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonded.add((u, v))
        bonded.add((v, u))
    thresh = 0.64  # 0.8² Å²
    score = 0
    for i in range(n):
        xi, yi = pos[i]
        for j in range(i + 1, n):
            if (i, j) not in bonded:
                dx, dy = xi - pos[j][0], yi - pos[j][1]
                if dx * dx + dy * dy < thresh:
                    score += 1
    return score


def _best_layout(romol, base_seed: int, n_tries: int = 6):
    """Try n_tries random atom orderings with CoordGen; apply the layout with fewest
    overlaps back to romol using the original atom indices (preserves residue map)."""
    import random

    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.rdchem import Conformer

    n = romol.GetNumAtoms()
    best_positions = None
    best_score = 10**9

    # Compute baseline CoordGen score so we never return something worse
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(romol)
    baseline_score = _overlap_score(romol)
    baseline_conf = romol.GetConformer()
    baseline_positions = {
        i: (baseline_conf.GetAtomPosition(i).x, baseline_conf.GetAtomPosition(i).y)
        for i in range(n)
    }
    best_score = baseline_score
    best_positions = baseline_positions

    for i in range(n_tries):
        rng = random.Random(base_seed + i)
        new_order = list(range(n))
        rng.shuffle(new_order)
        candidate = Chem.RenumberAtoms(romol, new_order)
        try:
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(candidate)
        except Exception:
            continue
        score = _overlap_score(candidate)
        if score < best_score:
            best_score = score
            conf = candidate.GetConformer()
            # new_order[new_idx] == old_idx → translate back to original indexing
            best_positions = {}
            for new_idx in range(n):
                p = conf.GetAtomPosition(new_idx)
                best_positions[new_order[new_idx]] = (p.x, p.y)

    if best_positions is None:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)
        return romol

    if romol.GetNumConformers() == 0:
        conf = Conformer(n)
        for idx, (x, y) in best_positions.items():
            conf.SetAtomPosition(idx, (x, y, 0.0))
        romol.AddConformer(conf, assignId=True)
    else:
        conf = romol.GetConformer()
        for idx, (x, y) in best_positions.items():
            conf.SetAtomPosition(idx, (x, y, 0.0))
    return romol


def _draw_mol(
    romol, width: int, height: int, used_slots: set | None = None, seed: int = 0
) -> str:
    import re as _re

    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    if seed == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)
    elif seed % 2 == 1:
        if _indigo_layout(romol) is None:
            romol = _best_layout(romol, base_seed=seed * 6)
    else:
        romol = _best_layout(romol, base_seed=seed * 6)
    rdDepictor.NormalizeDepiction(romol)
    rdDepictor.StraightenDepiction(romol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.padding = 0.12
    natoms = romol.GetNumAtoms()
    if natoms > 150:
        opts.minFontSize = 5
        opts.maxFontSize = 8
        opts.bondLineWidth = 0.8
        opts.additionalAtomLabelPadding = 0.0
    elif natoms > 80:
        opts.minFontSize = 7
        opts.maxFontSize = 10
        opts.bondLineWidth = 1.0
    drawer.DrawMolecule(romol)
    drawer.FinishDrawing()
    svg = _tag_atom_joins(drawer.GetDrawingText(), drawer, natoms)
    svg = _tag_stereo_annotations(svg, romol)

    if used_slots is not None:
        dummy_slot = {}
        for atom in romol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
                dummy_slot[atom.GetIdx()] = atom.GetIsotope()

        COLOR_USED = "#d9534f"
        COLOR_FREE = "#3dbe6c"
        OPACITY_USED = "0.35"

        for aidx, slot in dummy_slot.items():
            color = COLOR_USED if slot in used_slots else COLOR_FREE
            tag = f"atom-{aidx}"
            svg = _re.sub(
                rf"(<path\s+class='[^']*{tag}[^']*'[^>]*fill=')[^']+(')",
                rf"\g<1>{color}\2",
                svg,
            )
            if slot in used_slots:
                svg = _re.sub(
                    rf"(<path\s+class='[^']*{tag}[^']*'[^>]*)(/>)",
                    rf"\1 opacity='{OPACITY_USED}'\2",
                    svg,
                )
            svg = _re.sub(
                rf"(<path\s+class='bond-\d+\s+[^']*{tag}[^']*'[^>]*stroke:)#[0-9a-fA-F]{{6}}",
                rf"\g<1>{color}",
                svg,
            )
            if slot in used_slots:
                svg = _re.sub(
                    rf"(<path\s+class='bond-\d+\s+[^']*{tag}[^']*'[^>]*stroke-opacity:)\d[\d.]*",
                    rf"\g<1>{OPACITY_USED}",
                    svg,
                )
    return svg


def _mol_block(romol) -> str:
    from rdkit.Chem import MolToMolBlock

    return MolToMolBlock(romol)


def _canon(romol) -> str:
    from rdkit.Chem import MolToSmiles

    return MolToSmiles(romol, canonical=True)


def _count_defined_stereo(mol) -> int:
    """Count supplied atom and bond stereochemistry, including E/Z double bonds."""
    from rdkit import Chem

    return sum(
        1 for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ) + sum(
        bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
        for bond in mol.GetBonds()
    )


def _canon_flat(romol) -> str:
    """Canonical SMILES with all stereo stripped."""
    from rdkit import Chem

    rw = Chem.RWMol(romol)
    Chem.RemoveStereochemistry(rw)
    return Chem.MolToSmiles(rw, canonical=True)
