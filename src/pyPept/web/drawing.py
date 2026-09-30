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


def _draw_svg(drawer, molecule):
    """Draw a prepared copy; map added display hydrogens to their source parent."""
    import re

    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    source_index = "_cabiln_draw_source_atom"
    marked = Chem.Mol(molecule)
    for atom in marked.GetAtoms():
        atom.SetIntProp(source_index, atom.GetIdx())
    prepared = rdMolDraw2D.PrepareMolForDrawing(marked)
    drawer.drawOptions().prepareMolsBeforeDrawing = False
    drawer.DrawMolecule(prepared)
    drawer.FinishDrawing()
    svg = _tag_atom_joins(drawer.GetDrawingText(), drawer, prepared.GetNumAtoms())
    svg = _tag_stereo_annotations(svg, prepared)

    owners = {}
    for atom in prepared.GetAtoms():
        source = atom
        if (
            not atom.HasProp(source_index)
            and atom.GetAtomicNum() == 1 and atom.GetDegree() == 1
        ):
            source = atom.GetNeighbors()[0]
        if source.HasProp(source_index):
            owner = source.GetIntProp(source_index)
            if owner != atom.GetIdx():
                owners[atom.GetIdx()] = owner
    if not owners:
        return svg
    mapped_classes = {f"atom-{index}": f"atom-{owner}" for index, owner in owners.items()}

    def remap_classes(match):
        classes = (mapped_classes.get(name, name) for name in match[1].split())
        return "class='" + " ".join(dict.fromkeys(classes)) + "'"

    return re.sub(r"class='([^']*)'", remap_classes, svg)


def _indigo_layout(romol):
    """Accept Indigo coordinates only for the same mapped, stereochemical graph."""
    import math

    from rdkit import Chem

    try:
        # Plain SMILES cannot carry queries or every RDKit stereo/bond type.
        if any(
            atom.HasQuery() or atom.GetAtomicNum() == 0 or atom.GetChiralTag() not in (
                Chem.ChiralType.CHI_UNSPECIFIED, Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            )
            for atom in romol.GetAtoms()
        ) or any(
            bond.HasQuery() or bond.GetBondType() not in (
                Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                Chem.BondType.TRIPLE, Chem.BondType.AROMATIC,
            ) or bond.GetStereo() not in (
                Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY,
                Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS,
                Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
            )
            for bond in romol.GetBonds()
        ):
            return None
        from indigo import Indigo

        # A detached, mapped SMILES needs no preliminary coordinate calculation.
        # The original atoms, metadata, maps and stereo groups never leave RDKit.
        temporary = Chem.Mol(romol)
        for atom in temporary.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        expected = Chem.MolToSmiles(temporary)
        indigo = Indigo()
        im = indigo.loadMolecule(Chem.MolToSmiles(temporary, canonical=False))
        im.layout()
        result = Chem.MolFromMolBlock(im.molfile(), removeHs=False)
        if result is None or result.GetNumConformers() == 0:
            return None
        n = romol.GetNumAtoms()
        maps = [atom.GetAtomMapNum() for atom in result.GetAtoms()]
        if sorted(maps) != list(range(1, n + 1)) or Chem.MolToSmiles(result) != expected:
            return None
        source = result.GetConformer()
        positions = []
        for i in range(n):
            point = source.GetAtomPosition(i)
            positions.append((point.x, point.y, point.z))
        if source.Is3D() or any(
            not all(math.isfinite(value) for value in position) or position[2] != 0
            for position in positions
        ) or len(set(positions)) != n:
            return None
        conformer = Chem.Conformer(n)
        conformer.Set3D(False)
        for identifier, position in zip(maps, positions):
            conformer.SetAtomPosition(identifier - 1, position)
    except Exception:
        return None
    # Commit only after every external atom and coordinate passed validation.
    romol.RemoveAllConformers()
    romol.AddConformer(conformer, assignId=True)
    return romol


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
    # Compute baseline CoordGen score so we never return something worse
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(romol)
    baseline_score = _overlap_score(romol)
    if baseline_score == 0:
        return romol
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
            if score == 0:
                break

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
        if _indigo_layout(romol) is None:
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
    svg = _draw_svg(drawer, romol)

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
