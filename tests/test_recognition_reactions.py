"""Reaction proposals need source provenance and independent forward evidence."""

import pytest
from rdkit import Chem

from pyPept.interfaces.reaction_library import REACTIONS, run_bond_smirks
from pyPept.molecule import Molecule
from pyPept.recognition_reactions import recognition_variants
from pyPept.sequence import Sequence


def _assemble(notation):
    return Molecule(Sequence(notation, warning_sink=lambda _: None)).get_molecule(
        fmt="ROMol"
    )


def _smiles(molecule):
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def _forward(variant):
    """Reintroduce the declared attachment handles and run actual assembly."""
    assert len(variant.junctions) == 1
    junction = variant.junctions[0]
    molecule = Chem.RWMol(variant.molecule)
    for endpoint, isotope in ((junction.endpoint_a, 997), (junction.endpoint_b, 998)):
        anchor = molecule.GetAtomWithIdx(endpoint)
        if anchor.GetNumExplicitHs():
            anchor.SetNumExplicitHs(anchor.GetNumExplicitHs() - 1)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(isotope)
        index = molecule.AddAtom(dummy)
        molecule.AddBond(endpoint, index, Chem.BondType.SINGLE)
    Chem.SanitizeMol(molecule)
    return run_bond_smirks(
        molecule,
        molecule,
        997,
        998,
        REACTIONS[junction.reaction_id],
        intramolecular=True,
    )


@pytest.mark.parametrize(
    "product,precursor,reaction_id,virtual_count",
    [
        (
            "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am",
            "ac-Pra-A-A-AzK-am",
            "cuaac_1_4_triazole",
            0,
        ),
        (
            "ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am",
            "ac-S5-A-A-R8-am",
            "rcm_alkene",
            2,
        ),
    ],
)
@pytest.mark.parametrize("reverse_atom_order", [False, True])
def test_proposals_restore_precursor_and_reassemble_exact_product(
    product,
    precursor,
    reaction_id,
    virtual_count,
    reverse_atom_order,
):
    source = _assemble(product)
    if reverse_atom_order:
        source = Chem.RenumberAtoms(source, list(reversed(range(source.GetNumAtoms()))))
    for atom in source.GetAtoms():
        atom.SetIntProp("source_index_test", atom.GetIdx())
    original_smiles = _smiles(source)
    result = recognition_variants(source)
    assert not result.truncated
    assert not result.warnings
    assert _smiles(result.variants[0].molecule) == original_smiles
    assert result.variants[0].junctions == ()
    variant = next(
        v
        for v in result.variants
        if len(v.junctions) == 1 and v.junctions[0].reaction_id == reaction_id
    )
    assert _smiles(variant.molecule) == _smiles(_assemble(precursor))
    assert _smiles(_forward(variant)) == original_smiles
    assert (
        variant.atom_origins
        == tuple(range(source.GetNumAtoms())) + (None,) * virtual_count
    )
    junction = variant.junctions[0]
    assert len(junction.virtual_atoms) == virtual_count
    assert all(variant.atom_origins[i] is None for i in junction.virtual_atoms)
    assert {junction.endpoint_a, junction.endpoint_b} <= junction.source_atoms
    assert max(junction.source_atoms) < source.GetNumAtoms()
    assert junction.byproduct_smiles == ("C=C" if virtual_count else None)
    assert variant.changes[0].reaction_id == reaction_id
    assert variant.changes[0].bonds
    # Copies preserve original source identity even where bond/charge state changes.
    assert [
        a.GetIntProp("source_index_test")
        for a in variant.molecule.GetAtoms()
        if a.HasProp("source_index_test")
    ] == list(range(source.GetNumAtoms()))
    assert _smiles(source) == original_smiles
    assert source.GetNumAtoms() + virtual_count == variant.molecule.GetNumAtoms()


@pytest.mark.parametrize(
    "product,precursor,reaction_id,virtual_count",
    [
        (
            "ac-Pra.!1(4,4)-Pra.!2(4,4)-AzK.!1(4,4)-AzK.!2(4,4)-am",
            "ac-Pra-Pra-AzK-AzK-am",
            "cuaac_1_4_triazole",
            0,
        ),
        (
            "ac-S5.!1(4,4)-S5.!2(4,4)-A-R8.!1(4,4)-R8.!2(4,4)-am",
            "ac-S5-S5-A-R8-R8-am",
            "rcm_alkene",
            4,
        ),
    ],
)
def test_multiple_reaction_sites_are_not_truncated_to_first_match(
    product,
    precursor,
    reaction_id,
    virtual_count,
):
    source = _assemble(product)
    result = recognition_variants(source)
    variant = next(
        v
        for v in result.variants
        if len(v.junctions) == 2
        and all(j.reaction_id == reaction_id for j in v.junctions)
    )
    assert not result.truncated
    assert _smiles(variant.molecule) == _smiles(_assemble(precursor))
    assert sum(origin is None for origin in variant.atom_origins) == virtual_count
    assert len({(j.endpoint_a, j.endpoint_b) for j in variant.junctions}) == 2


def test_spaac_internal_alkyne_does_not_gain_a_terminal_hydrogen():
    # Independent assembly models, without monomer names or the peptide parser.
    cyclooctyne = Chem.MolFromSmiles("C1CC([400*])C#CCCC1")
    azide = Chem.MolFromSmiles("[401*]CN=[N+]=[N-]")
    source = run_bond_smirks(
        cyclooctyne,
        azide,
        400,
        401,
        REACTIONS["spaac_triazole"],
        intramolecular=False,
    )
    result = recognition_variants(source)
    variant = next(
        v
        for v in result.variants
        if any(j.reaction_id == "spaac_triazole" for j in v.junctions)
    )
    expected = Chem.MolFromSmiles("C1CCC#CCCC1.CN=[N+]=[N-]")
    assert _smiles(variant.molecule) == _smiles(expected)
    assert _smiles(_forward(variant)) == _smiles(source)
    assert variant.atom_origins == tuple(range(source.GetNumAtoms()))
    assert not any(
        j.reaction_id == "cuaac_1_4_triazole"
        for v in result.variants
        for j in v.junctions
    )


@pytest.mark.parametrize("smiles", ["N[C@@H](C)C(=O)NCC(=O)O", "Cc1nncn1C"])
def test_unrelated_molecules_are_not_normalized(smiles):
    source = Chem.MolFromSmiles(smiles)
    assert source is not None
    result = recognition_variants(source)
    assert len(result.variants) == 1
    assert not result.truncated
    assert _smiles(result.variants[0].molecule) == _smiles(source)
    assert result.variants[0].molecule is not source


@pytest.mark.parametrize("limit", [{"max_proposals": 1}, {"max_matches": 1}])
def test_budget_exhaustion_retains_original_with_explicit_status(limit):
    source = _assemble("ac-Pra.!1(4,4)-Pra.!2(4,4)-AzK.!1(4,4)-AzK.!2(4,4)-am")
    result = recognition_variants(source, **limit)
    assert result.truncated
    assert result.warnings
    assert _smiles(result.variants[0].molecule) == _smiles(source)


def test_declared_metathesis_byproduct_is_required(monkeypatch):
    from pyPept import recognition_reactions

    entry = dict(REACTIONS["rcm_alkene"])
    entry["steps"] = ["[4*][C:1]=[CH2:2].[4*][C:3]=[CH2:4] >> [C:1]=[C:3]"]
    monkeypatch.setattr(
        recognition_reactions, "REACTIONS", {"altered_metathesis": entry}
    )
    source = _assemble("ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am")
    assert len(recognition_variants(source).variants) == 1


def test_declared_types_and_product_rule_work_without_reaction_name_whitelist(
    monkeypatch,
):
    from pyPept import recognition_reactions

    monkeypatch.setattr(
        recognition_reactions, "REACTIONS", {"renamed_rule": REACTIONS["rcm_alkene"]}
    )
    source = _assemble("ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am")
    variants = recognition_variants(source).variants
    assert variants[1].junctions[0].reaction_id == "renamed_rule"
    assert _smiles(variants[1].molecule) == _smiles(_assemble("ac-S5-A-A-R8-am"))
