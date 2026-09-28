"""Backbone preference must preserve chemically defined monomer boundaries."""

import pytest
from rdkit import Chem

from pyPept import monomer_store
from pyPept.molecule import Molecule
from pyPept.recognition import recognize
from pyPept.recognition_notation import emit_interpretation
from pyPept.sequence import Sequence
from pyPept.smiles import convert_smiles

# Literal input atom groups and junctions are annotated from the supplied SMILES,
# independently of CABILN assembly. Every listed bond is directed R2 -> R1;
# the flags distinguish whole backbone monomers from terminal caps.
CASES = [
    pytest.param(
        "CNCC(=O)NCC(=O)O",
        (range(5), range(5, 10)),
        ((3, 5),),
        (True, True),
        1,
        id="sarcosine-glycine",
    ),
    pytest.param(
        "CN[C@@H](C)C(=O)NCC(=O)O",
        (range(6), range(6, 11)),
        ((4, 6),),
        (True, True),
        1,
        id="N-methyl-alanine-glycine",
    ),
    pytest.param(
        "NCC(=O)N(C)CC(=O)NCC(=O)O",
        (range(4), range(4, 9), range(9, 14)),
        ((2, 4), (7, 9)),
        (True, True, True),
        2,
        id="internal-N-methyl",
    ),
    pytest.param(
        "NCC(=O)NCC(=O)OC",
        (range(4), range(4, 8), range(8, 10)),
        ((2, 4), (6, 8)),
        (True, True, False),
        1,
        id="O-methyl-ester-cap",
    ),
    pytest.param(
        "N[C@@H](Cc1ccc(OC)cc1)C(=O)NCC(=O)O",
        (range(13), range(13, 18)),
        ((11, 13),),
        (True, True),
        1,
        id="O-methyl-sidechain",
    ),
    pytest.param(
        "CC(=O)N[C@@H](C)C(=O)NCC(=O)N",
        (range(3), range(3, 8), range(8, 12), (12,)),
        ((1, 3), (6, 8), (10, 12)),
        (False, True, True, False),
        1,
        id="acetyl-and-amide-caps",
    ),
    pytest.param(
        "NCC(=O)NCCNCC(=O)O",
        (range(4), range(4, 7), range(7, 12)),
        ((2, 4), (6, 7)),
        (True, True, True),
        2,
        id="reduced-amide-backbone",
    ),
]


def assert_assembly_and_ownership(source, notation, assignments):
    assembled = Molecule(
        Sequence(notation, warning_sink=lambda _: None), depiction=None
    )
    rebuilt = assembled.mol
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(source)
    assert sorted(atom for item in assignments for atom in item.source_atoms) == list(
        range(source.GetNumAtoms())
    )
    assert all(item.recognized for item in assignments)
    owned = assembled.get_residue_atom_map()
    matches = source.GetSubstructMatches(
        rebuilt, useChirality=True, uniquify=False, maxMatches=1000
    )
    assert 0 < len(matches) < 1000
    assert any(
        all(
            {mapping[atom] for atom in owned[item.residue_index]}
            == set(item.source_atoms)
            for item in assignments
        )
        for mapping in matches
    ), "The emitted occurrences must own the atoms reported by recognition"


@pytest.mark.parametrize("admitted", [False, True], ids=["proposal", "admitted"])
@pytest.mark.parametrize("smiles,groups,bonds,backbones,backbone_count", CASES)
def test_preference_retains_whole_monomers_and_real_junctions(
    smiles, groups, bonds, backbones, backbone_count, admitted
):
    source = Chem.MolFromSmiles(smiles)
    options = {"accept_cover": lambda _: True} if admitted else {}
    cover = recognize(source, **options).covers[0]
    by_atoms = {candidate.atoms: candidate for candidate in cover}
    assert set(by_atoms) == {frozenset(group) for group in groups}
    for group, backbone in zip(groups, backbones):
        assert by_atoms[frozenset(group)].has_backbone is backbone
    expected_ports = {(a, b, 2) for a, b in bonds} | {(b, a, 1) for a, b in bonds}
    assert {port for candidate in cover for port in candidate.ports} == expected_ports
    for candidate in cover:
        for inside, outside, slot in candidate.ports:
            assert inside in candidate.atoms and outside not in candidate.atoms
            assert (inside, slot) in candidate.anchors
            assert source.GetBondBetweenAtoms(inside, outside) is not None
    # Carbonyl-free reduced exits remain legitimate backbone sites.
    if smiles == "NCC(=O)NCCNCC(=O)O":
        assert dict(by_atoms[frozenset(range(4, 7))].attachment_types)[2] == (
            "backbone_c_red"
        )
    interpretation = emit_interpretation(source, cover)
    assert interpretation.recognized_atoms == source.GetNumAtoms()
    assert interpretation.backbone_connections == backbone_count
    assert_assembly_and_ownership(
        source, interpretation.cabiln, interpretation.assignments
    )
    result = convert_smiles(smiles)
    assert result.recognition_status == "complete"
    assert {frozenset(item.source_atoms) for item in result.assignments} == set(
        by_atoms
    )
    assert_assembly_and_ownership(source, result.cabiln, result.assignments)


@pytest.fixture
def small_library(tmp_path, monkeypatch):
    _, bundled = monomer_store._load_sdf()

    def write(renames):
        path = tmp_path / "monomers.sdf"
        with Chem.SDWriter(str(path)) as writer:
            for original, new in renames.items():
                molecule = Chem.Mol(bundled[original])
                molecule.SetProp("symbol", new)
                molecule.SetProp("m_abbr", new)
                writer.write(molecule)
        monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))

    return write


def test_preference_uses_structure_after_arbitrary_library_renaming(small_library):
    small_library({"meG": "same7New", "G": "DifferentGly", "Me_": "NewCap"})
    source = Chem.MolFromSmiles("CNCC(=O)NCC(=O)O")
    result = convert_smiles(Chem.MolToSmiles(source, canonical=False))
    assert result.recognition_status == "complete"
    assert {frozenset(item.source_atoms) for item in result.assignments} == {
        frozenset(range(5)),
        frozenset(range(5, 10)),
    }
    assert {item.symbol for item in result.assignments} == {"same7New", "DifferentGly"}
    assert_assembly_and_ownership(source, result.cabiln, result.assignments)


def test_methyl_cap_remains_available_when_whole_monomer_is_absent(small_library):
    small_library({"G": "G", "Me_": "Me_"})
    source = Chem.MolFromSmiles("CNCC(=O)NCC(=O)O")
    cover = recognize(source).covers[0]
    assert {candidate.atoms for candidate in cover} == {
        frozenset((0,)),
        frozenset(range(1, 5)),
        frozenset(range(5, 10)),
    }
    cap = next(candidate for candidate in cover if candidate.atoms == frozenset((0,)))
    assert cap.ports == ((0, 1, 2),)
    assert dict(cap.attachment_types)[2] == "alkyl_halide_c"
    assert not cap.has_backbone
    interpretation = emit_interpretation(source, cover)
    assert interpretation.backbone_connections == 1
    assert interpretation.recognized_atoms == source.GetNumAtoms()
    assert_assembly_and_ownership(
        source, interpretation.cabiln, interpretation.assignments
    )
