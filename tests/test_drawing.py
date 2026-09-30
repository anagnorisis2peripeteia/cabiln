"""Depiction coordinates and SVG identity used by highlighting and selection."""

import math
import re
import sys
from types import SimpleNamespace
from unittest.mock import Mock
from xml.etree import ElementTree as ET

import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from pyPept.inputs import read_input
from pyPept.web import drawing
from pyPept.web.drawing import _draw_mol, _tag_atom_joins, _tag_stereo_annotations


def _without_coordinates(molecule):
    copied = Chem.Mol(molecule)
    copied.RemoveAllConformers()
    return copied.ToBinary(Chem.PropertyPickleOptions.AllProps)


@pytest.mark.parametrize("source,kind", [
    ("C-S-K.[G(4,2)]-am", "cabiln"),
    ("[NH3+][13CH2]C(=O)O.C/C=C/C.C/C=C\\C", "smiles"),
    ("N[C@@H](C)C(=O)N[C@@H](C)C(=O)O |a:1,6|", "smiles"),
])
def test_indigo_layout_preserves_chemistry_and_ownership_without_bootstrap(
    monkeypatch, source, kind
):
    pytest.importorskip("indigo")
    molecule = read_input(source, input_format=kind).assemble(depiction=None)
    molecule.RemoveAllConformers()
    assert molecule.GetNumConformers() == 0
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 31)
    molecule.SetProp("source", source)
    molecule.GetBondWithIdx(0).SetProp("bondNote", "preserve this")
    expected = _without_coordinates(molecule)
    expected_smiles = Chem.MolToSmiles(Chem.Mol(molecule))
    bootstrap = Mock(side_effect=AssertionError("Discarded preliminary layout"))
    monkeypatch.setattr(rdDepictor, "Compute2DCoords", bootstrap)

    assert drawing._indigo_layout(molecule) is molecule
    bootstrap.assert_not_called()
    assert _without_coordinates(molecule) == expected
    conformer = molecule.GetConformer()
    assert not conformer.Is3D()
    positions = [tuple(conformer.GetAtomPosition(i)) for i in range(molecule.GetNumAtoms())]
    assert all(math.isfinite(value) for position in positions for value in position)
    assert all(position[2] == 0 for position in positions)
    assert len(set(positions)) == len(positions)
    exported = Chem.MolFromMolBlock(drawing._mol_block(molecule), removeHs=False)
    assert Chem.MolToSmiles(exported) == expected_smiles


def _external_layout(monkeypatch, returned):
    """Replace the external result while retaining the real acceptance checks."""
    native = Mock()
    native.molfile.return_value = "external layout"
    session = Mock()
    session.loadMolecule.return_value = native
    factory = Mock(return_value=session)
    monkeypatch.setitem(sys.modules, "indigo", SimpleNamespace(Indigo=factory))
    monkeypatch.setattr(Chem, "MolFromMolBlock", Mock(return_value=returned))
    return factory


@pytest.mark.parametrize("fault", [
    "missing-map", "duplicate-map", "extra-atom", "bond", "atom-stereo",
    "bond-stereo", "nonfinite", "3d", "collapsed", "unreadable", "no-conformer",
])
def test_invalid_external_layout_leaves_source_intact_and_uses_coordgen(monkeypatch, fault):
    molecule = Chem.MolFromSmiles("N[C@@H](C)C(=O)O.C/C=C/C")
    rdDepictor.Compute2DCoords(molecule, forceRDKit=True)
    before = molecule.ToBinary(Chem.PropertyPickleOptions.AllProps)
    expected_smiles = Chem.MolToSmiles(Chem.Mol(molecule))
    returned = Chem.RWMol(molecule)
    for atom in returned.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    if fault == "missing-map":
        returned.GetAtomWithIdx(0).SetAtomMapNum(0)
    elif fault == "duplicate-map":
        returned.GetAtomWithIdx(0).SetAtomMapNum(2)
    elif fault == "extra-atom":
        returned.AddAtom(Chem.Atom(6))
    elif fault == "bond":
        returned.GetBondBetweenAtoms(6, 7).SetBondType(Chem.BondType.DOUBLE)
    elif fault == "atom-stereo":
        returned.GetAtomWithIdx(1).InvertChirality()
    elif fault == "bond-stereo":
        returned.GetBondBetweenAtoms(7, 8).SetStereo(Chem.BondStereo.STEREOZ)
        Chem.SetDoubleBondNeighborDirections(returned)
    elif fault == "nonfinite":
        returned.GetConformer().SetAtomPosition(0, (math.nan, 0, 0))
    elif fault == "3d":
        returned.GetConformer().Set3D(True)
        returned.GetConformer().SetAtomPosition(0, (0, 0, 1))
    elif fault == "collapsed":
        for i in range(returned.GetNumAtoms()):
            returned.GetConformer().SetAtomPosition(i, (0, 0, 0))
    elif fault == "unreadable":
        returned = None
    elif fault == "no-conformer":
        returned.RemoveAllConformers()
    _external_layout(monkeypatch, returned)

    assert drawing._indigo_layout(molecule) is None
    assert molecule.ToBinary(Chem.PropertyPickleOptions.AllProps) == before
    fallback = Mock(wraps=rdDepictor.Compute2DCoords)
    monkeypatch.setattr(rdDepictor, "Compute2DCoords", fallback)
    assert "<svg" in _draw_mol(molecule, 400, 300)
    assert fallback.call_count >= 1
    assert Chem.MolToSmiles(molecule) == expected_smiles


def test_external_atom_reordering_copies_coordinates_by_map_id(monkeypatch):
    molecule = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    returned = Chem.Mol(molecule)
    for atom in returned.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    rdDepictor.Compute2DCoords(returned, forceRDKit=True)
    positions = [tuple(p) for p in returned.GetConformer().GetPositions()]
    returned = Chem.RenumberAtoms(returned, list(reversed(range(returned.GetNumAtoms()))))
    factory = _external_layout(monkeypatch, returned)
    before = _without_coordinates(molecule)

    assert drawing._indigo_layout(molecule) is molecule
    assert [tuple(p) for p in molecule.GetConformer().GetPositions()] == positions
    assert _without_coordinates(molecule) == before
    transport = factory.return_value.loadMolecule.call_args.args[0]
    mapped = Chem.MolFromSmiles(transport)
    assert [atom.GetAtomMapNum() for atom in mapped.GetAtoms()] == list(range(1, 7))


@pytest.mark.parametrize("availability", ["missing", "constructor-failure", "layout-failure"])
def test_unavailable_indigo_falls_back_to_coordgen(monkeypatch, availability):
    factory = _external_layout(monkeypatch, None)
    if availability == "missing":
        monkeypatch.setitem(sys.modules, "indigo", None)
    elif availability == "constructor-failure":
        factory.side_effect = RuntimeError("native library unavailable")
    else:
        factory.return_value.loadMolecule.return_value.layout.side_effect = RuntimeError("layout failed")
    molecule = Chem.MolFromSmiles("NCC(=O)O")
    fallback = Mock(wraps=rdDepictor.Compute2DCoords)
    monkeypatch.setattr(rdDepictor, "Compute2DCoords", fallback)
    assert "<svg" in _draw_mol(molecule, 400, 300)
    assert fallback.call_count >= 1


@pytest.mark.parametrize("kind", ["dummy", "query", "dative", "square-planar", "atropisomer"])
def test_unsupported_transport_skips_indigo(monkeypatch, kind):
    molecules = {
        "dummy": Chem.MolFromSmiles("[1*]NCC(=O)O"),
        "query": Chem.MolFromSmarts("[#6]~[#7]"),
        "dative": Chem.MolFromSmiles("N->[Cu+2]"),
        "square-planar": Chem.MolFromSmiles("N[Pt@SP1](N)(Cl)Cl"),
        "atropisomer": Chem.MolFromSmiles("Cc1cccc(F)c1-c1c(Cl)cccc1Br"),
    }
    molecule = molecules[kind]
    if kind == "atropisomer":
        axis = next(bond for bond in molecule.GetBonds()
                    if not bond.IsInRing() and bond.GetBeginAtom().GetIsAromatic()
                    and bond.GetEndAtom().GetIsAromatic())
        axis.SetStereo(Chem.BondStereo.STEREOATROPCW)
    factory = _external_layout(monkeypatch, None)
    before = molecule.ToBinary(Chem.PropertyPickleOptions.AllProps)
    assert drawing._indigo_layout(molecule) is None
    factory.assert_not_called()
    assert molecule.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("scores,tries,chosen", [
    ([0], 6, 1), ([3, 0], 6, 2), ([3, 2, 2], 2, 2),
], ids=["zero-baseline", "zero-candidate", "keep-earlier-tie"])
def test_coordgen_search_stops_at_zero_and_retains_strict_ties(monkeypatch, scores, tries, chosen):
    molecule = Chem.MolFromSmiles("CCCC")
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    layouts = []

    def coordinates(mol):
        layouts.append(mol)
        conformer = Chem.Conformer(mol.GetNumAtoms())
        conformer.Set3D(False)
        for atom in mol.GetAtoms():
            conformer.SetAtomPosition(atom.GetIdx(), (10 * len(layouts) + atom.GetAtomMapNum(), 0, 0))
        mol.RemoveAllConformers()
        mol.AddConformer(conformer)

    monkeypatch.setattr(rdDepictor, "Compute2DCoords", coordinates)
    monkeypatch.setattr(drawing, "_overlap_score", Mock(side_effect=scores))
    assert drawing._best_layout(molecule, base_seed=12, n_tries=tries) is molecule
    assert len(layouts) == len(scores)
    assert [p[0] for p in molecule.GetConformer().GetPositions()] == [
        chosen * 10 + i + 1 for i in range(4)
    ]


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


@pytest.mark.parametrize("explicit_hydrogen", [False, True], ids=["drawing-hydrogen", "source-hydrogen"])
def test_drawing_hydrogen_uses_parent_ownership_without_changing_source(explicit_hydrogen):
    source = "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3"
    molecule = read_input(source).assemble(depiction=None)
    if explicit_hydrogen:
        molecule = Chem.AddHs(molecule, onlyOnAtoms=[19])
    source_atoms = molecule.GetNumAtoms()
    source_bonds = molecule.GetNumBonds()
    expected_smiles = Chem.MolToSmiles(Chem.Mol(molecule))
    expected_properties = _without_coordinates(molecule)

    svg = _draw_mol(molecule, 960, 680)
    assert molecule.GetNumAtoms() == source_atoms
    assert molecule.GetNumBonds() == source_bonds
    assert _without_coordinates(molecule) == expected_properties
    exported = Chem.MolFromMolBlock(drawing._mol_block(molecule), removeHs=False)
    assert Chem.MolToSmiles(exported) == expected_smiles

    prepared = rdMolDraw2D.PrepareMolForDrawing(molecule)
    hydrogen = prepared.GetAtomWithIdx(51)
    assert hydrogen.GetAtomicNum() == 1
    assert [atom.GetIdx() for atom in hydrogen.GetNeighbors()] == [19]
    assert prepared.GetNumAtoms() == 52
    expected_owner = 51 if explicit_hydrogen else 19

    # Native automatic preparation supplies an independent geometry/class oracle.
    drawer = rdMolDraw2D.MolDraw2DSVG(960, 680)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().padding = 0.12
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    paths = ET.fromstring(svg).findall("{http://www.w3.org/2000/svg}path")
    native_paths = ET.fromstring(drawer.GetDrawingText()).findall("{http://www.w3.org/2000/svg}path")
    assert len(paths) == len(native_paths)
    labels, bonds = 0, 0
    for path, native in zip(paths, native_paths):
        assert {k: v for k, v in path.attrib.items() if k != "class"} == {
            k: v for k, v in native.attrib.items() if k != "class"
        }
        classes = native.get("class", "").split()
        if "atom-51" not in classes:
            continue
        owners = {int(value[5:]) for value in path.get("class", "").split()
                  if value.startswith("atom-")}
        if any(value.startswith("bond-") for value in classes):
            bonds += 1
            assert owners == {19, expected_owner}
        else:
            labels += 1
            assert owners == {expected_owner}
    assert labels and bonds


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
