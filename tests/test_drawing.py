"""SVG atom identity used by residue highlighting and pointer selection."""

import re
from xml.etree import ElementTree as ET

import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from pyPept.inputs import read_input
from pyPept.web.drawing import _draw_mol, _tag_atom_joins, _tag_stereo_annotations


@pytest.mark.parametrize("source,seed", [("G", 0), ("C-S-K.[G(4,2)]-am", 2)])
def test_atom_labels_and_bond_joins_keep_their_atom_identity(source, seed):
    parsed = read_input(source)
    molecule = parsed.assemble()
    paths = ET.fromstring(_draw_mol(molecule, 960, 680, seed=seed)).findall(
        "{http://www.w3.org/2000/svg}path"
    )
    owners = {
        index for indices in parsed.assembly.get_residue_atom_map().values()
        for index in indices
    }
    assert owners == set(range(molecule.GetNumAtoms()))

    joins = [
        path for path in paths
        if "stroke-miterlimit:10" in path.get("style", "")
        and path.get("style", "").startswith("fill:none;")
    ]
    assert joins
    for join in joins:
        classes = join.get("class", "").split()
        assert len(classes) == 1 and classes[0].startswith("atom-"), join.attrib
        # A join's middle point must be shared by two different bonds at its atom.
        vertex = re.findall(r"-?\d+\.\d+,-?\d+\.\d+", join.get("d"))[1]
        incident = {
            path.get("class").split()[0]
            for path in paths
            if path.get("class", "").startswith("bond-")
            and classes[0] in path.get("class").split()
            and vertex in re.findall(r"-?\d+\.\d+,-?\d+\.\d+", path.get("d"))
        }
        assert len(incident) >= 2, join.attrib

    for atom in molecule.GetAtoms():
        if atom.GetSymbol() not in ("N", "O", "S"):
            continue
        # NH2 consists of three glyphs; the implicit Hs belong to the same atom.
        glyphs = [
            path for path in paths
            if path.get("fill") and path.get("class") == f"atom-{atom.GetIdx()}"
        ]
        hydrogens = atom.GetTotalNumHs()
        assert len(glyphs) == 1 + bool(hydrogens) + (hydrogens > 1)


def test_coincident_atom_joins_do_not_claim_an_arbitrary_owner():
    molecule = Chem.MolFromSmiles("CCC.CCC")
    conformer = Chem.Conformer(6)
    for index, (x, y) in enumerate([(0, 0), (1, 1), (2, 0)] * 2):
        conformer.SetAtomPosition(index, (x, y, 0))
    molecule.AddConformer(conformer)
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    assert "<path d='M" in svg
    assert _tag_atom_joins(svg, drawer, molecule.GetNumAtoms()) == svg


@pytest.mark.parametrize("no_freetype", [False, True], ids=["paths", "text"])
def test_stereo_glyph_owners_match_independently_drawn_single_notes(no_freetype):
    molecule = Chem.MolFromSmiles(
        "N[C@@H](C)C(=O)NCCCCCC[C@H](F)Cl.C/C=C/C.C/C=C\\C"
    )
    molecule.GetAtomWithIdx(0).SetProp("atomNote", "custom")
    molecule.GetBondWithIdx(0).SetProp("bondNote", "extra")
    rdDepictor.Compute2DCoords(molecule)
    source = Chem.MolToCXSmiles(molecule)
    expected = [("atom", 1, {1}), ("atom", 12, {12}),
                ("bond", 15, {16, 17}), ("bond", 18, {20, 21})]

    def glyphs(mol):
        drawer = rdMolDraw2D.MolDraw2DSVG(960, 680, noFreetype=no_freetype)
        options = drawer.drawOptions()
        options.addStereoAnnotation = True
        # Hold geometry and scale fixed while removing unrelated annotations.
        options.fixedBondLength = 24
        options.fixedFontSize = 14
        options.drawingExtentsInclude = (
            rdMolDraw2D.DrawElement.ALL ^ rdMolDraw2D.DrawElement.ANNOTATIONS
        )
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        original = drawer.GetDrawingText()
        tagged = _tag_stereo_annotations(original, mol)
        assert re.sub(r"class='CIP_Code[^']*'", "class='CIP_Code'", tagged) == original
        unexpected = original.replace("class='CIP_Code'", "class='note'", 1)
        assert _tag_stereo_annotations(unexpected, mol) == unexpected
        root = ET.fromstring(tagged)
        return {
            (element.tag, tuple(sorted((k, v) for k, v in element.attrib.items()
                                      if k != "class")), element.text):
            set(element.get("class", "").split())
            for element in root if "CIP_Code" in element.get("class", "").split()
        }

    complete = glyphs(molecule)
    assert len(complete) == 12
    for kind, index, atoms in expected:
        isolated = Chem.Mol(molecule)
        for atom in isolated.GetAtoms():
            if (kind, index) != ("atom", atom.GetIdx()) and atom.HasProp("_CIPCode"):
                atom.ClearProp("_CIPCode")
        for bond in isolated.GetBonds():
            if (kind, index) != ("bond", bond.GetIdx()):
                bond.SetStereo(Chem.BondStereo.STEREONONE)
                if bond.HasProp("_CIPCode"):
                    bond.ClearProp("_CIPCode")
        single = glyphs(isolated)
        assert len(single) == 3
        assert single.keys() <= complete.keys()
        for signature in single:
            assert {value for value in complete[signature] if value.startswith("atom-")} == {
                f"atom-{atom}" for atom in atoms
            }
    assert Chem.MolToCXSmiles(molecule) == source
    assert molecule.GetAtomWithIdx(0).GetProp("atomNote") == "custom"
    assert molecule.GetBondWithIdx(0).GetProp("bondNote") == "extra"
