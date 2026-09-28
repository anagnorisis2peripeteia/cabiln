"""Acceptance checks for residue recovery, independent of molecule round-tripping.

Expected partitions come from separately supplied peptide constructions. A
coarse whole-component escape cannot satisfy the monomer/template assertions,
even if its assembly recovers every atom of the input molecule.
"""

import random
import shutil
from collections import Counter

import pytest
from rdkit import Chem

from pyPept.interfaces.cli_monomer import register_monomer
from pyPept.molecule import Molecule
from pyPept.monomer_store import _load_sdf, library_path
from pyPept.recognition import RecognitionBudgets
from pyPept.sequence import Sequence, _rgroup_atom_idx
from pyPept.smiles import convert_smiles

TWO_CHAINS = "C.!1(4,4)-G%C.!1(4,4)-G"
MIXED_BACKBONE = "ac-A-redG2-G-am"
BRANCH_AND_SECOND_CHAIN = "ac-K.[G(4,2).ac(1,2)]-C.!1(4,4)-am%ac-C.!1-G-am"
LIBRARY_ORDER_CASE = "ac-K.!1(4,2)-A-am%ac-G-G.!1(2,4)"
UNKNOWN = "<[1*]N[C@@H](CCC(F)(F)F)C([2*])=O>"
UNKNOWN_CONTEXT = f"A-G-{UNKNOWN}-K-V"


def _assembled(notation):
    sequence = Sequence(notation, warning_sink=lambda _message: None)
    assembly = Molecule(sequence)
    molecule = Chem.Mol(assembly.get_molecule("ROMol"))
    ownership = assembly.get_residue_atom_map()
    all_indices = [atom for atoms in ownership.values() for atom in atoms]
    assert len(all_indices) == len(set(all_indices)), "Atoms have multiple owners"
    assert set(all_indices) == set(range(molecule.GetNumAtoms())), "Unowned atoms"
    for residue, atoms in ownership.items():
        for atom in atoms:
            molecule.GetAtomWithIdx(atom).SetIntProp("_acceptance_owner", residue)
    molecule = Chem.RemoveHs(molecule)
    partition = [set() for _ in sequence.s_monomers]
    for atom in molecule.GetAtoms():
        assert atom.HasProp("_acceptance_owner"), "Source atom ownership was lost"
        partition[atom.GetIntProp("_acceptance_owner")].add(atom.GetIdx())
    assert all(partition), "An emitted residue owns no product atoms"
    return sequence, molecule, partition


def _template(monomer):
    # Keep actual numbered attachments and stereochemistry. Equivalent symbols
    # for the same bonded molecular template need not have identical names.
    return Chem.MolToSmiles(monomer["m_romol"], isomericSmiles=True)


def _synthetic_indices(sequence):
    synthetic = set(sequence._synthetic_aa_smiles) | set(sequence._synthetic_caps)
    return {
        i
        for i, monomer in enumerate(sequence.s_monomers)
        if monomer["m_abbr"] in synthetic
    }


def _partition_key(partition):
    return tuple(sorted(tuple(sorted(atoms)) for atoms in partition))


def _assert_exact_partition(expected, actual):
    _, source, source_partition = expected
    _, rebuilt, rebuilt_partition = actual
    assert Chem.MolToSmiles(rebuilt) == Chem.MolToSmiles(
        source
    ), "Molecular identity changed"
    mappings = source.GetSubstructMatches(
        rebuilt, useChirality=True, uniquify=False, maxMatches=1024
    )
    assert 0 < len(mappings) < 1024, "Incomplete molecular-isomorphism comparison"
    wanted = _partition_key(source_partition)
    assert any(
        _partition_key(
            [{mapping[atom] for atom in atoms} for atoms in rebuilt_partition]
        )
        == wanted
        for mapping in mappings
    ), "The product matches, but original residue atom boundaries were lost"


def _input_smiles(molecule, seed=None):
    if seed is None:
        return Chem.MolToSmiles(molecule, isomericSmiles=True)
    order = list(range(molecule.GetNumAtoms()))
    random.Random(seed).shuffle(order)
    renamed = Chem.RenumberAtoms(molecule, order)
    value = Chem.MolToRandomSmilesVect(
        renamed, 1, randomSeed=seed, isomericSmiles=True
    )[0]
    assert Chem.MolToSmiles(Chem.MolFromSmiles(value)) == Chem.MolToSmiles(molecule)
    return value


def _assert_registered_recovery(expected, result):
    actual = _assembled(result.cabiln)
    original_sequence, _, _ = expected
    recovered_sequence, _, _ = actual
    assert not _synthetic_indices(recovered_sequence), result.cabiln
    assert len(recovered_sequence.s_monomers) == len(
        original_sequence.s_monomers
    ), result.cabiln
    assert Counter(map(_template, recovered_sequence.s_monomers)) == Counter(
        map(_template, original_sequence.s_monomers)
    ), result.cabiln
    _assert_exact_partition(expected, actual)


@pytest.mark.parametrize(
    "notation",
    [
        pytest.param("A-G-A-G-A", id="repeated-residues"),
        pytest.param("ac-bAla-Aib-meNle-am", id="noncanonical-residues"),
        pytest.param(TWO_CHAINS, id="two-Cys-G-chains-disulfide"),
        pytest.param(MIXED_BACKBONE, id="amide-and-reduced-backbone"),
        pytest.param(BRANCH_AND_SECOND_CHAIN, id="branch-and-second-chain"),
        pytest.param(LIBRARY_ORDER_CASE, id="isopeptide-branch"),
    ],
)
def test_registered_peptides_keep_their_residue_partition(notation):
    expected = _assembled(notation)
    result = convert_smiles(_input_smiles(expected[1]))
    _assert_registered_recovery(expected, result)


@pytest.mark.parametrize("seed", [17, 41, 73])
@pytest.mark.parametrize("notation", ["A-G-A-G-A", TWO_CHAINS, "!1-A-G-K-G-!1"])
def test_registered_partition_is_independent_of_input_atom_order(notation, seed):
    expected = _assembled(notation)
    result = convert_smiles(_input_smiles(expected[1], seed))
    _assert_registered_recovery(expected, result)


@pytest.mark.parametrize("seed", [None, 17, 41, 73])
def test_unknown_internal_residue_keeps_known_neighbors_and_backbone_sites(seed):
    expected = _assembled(UNKNOWN_CONTEXT)
    result = convert_smiles(_input_smiles(expected[1], seed))
    actual = _assembled(result.cabiln)
    sequence = actual[0]
    synthetic = _synthetic_indices(sequence)
    assert len(sequence.s_monomers) == 5, result.cabiln
    assert len(synthetic) == 1, result.cabiln
    known_templates = Counter(
        _template(m) for i, m in enumerate(sequence.s_monomers) if i not in synthetic
    )
    assert known_templates == Counter(
        _template(m) for i, m in enumerate(expected[0].s_monomers) if i != 2
    ), result.cabiln
    local = sequence.s_monomers[next(iter(synthetic))]
    assert _template(local) == _template(
        expected[0].s_monomers[2]
    ), "The local unknown must retain its actual N/C backbone attachment sites"
    for slot, element in ((1, "N"), (2, "C")):
        dummy = _rgroup_atom_idx(local["m_romol"], slot)
        assert dummy is not None
        neighbors = local["m_romol"].GetAtomWithIdx(dummy).GetNeighbors()
        assert len(neighbors) == 1 and neighbors[0].GetSymbol() == element
    _assert_exact_partition(expected, actual)
    assert not any("coarse" in message.lower() for message in result.warnings)


def _copy_library(directory, seed=None):
    original = library_path()
    records = list(_load_sdf()[0])
    if seed is not None:
        random.Random(seed).shuffle(records)
    path = directory / "monomers.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for record in records:
            writer.write(record)
    companion = original.with_suffix(".csv")
    if companion.is_file():
        shutil.copyfile(companion, path.with_suffix(".csv"))
    return path


@pytest.mark.parametrize("seed", [17, 41, 73])
def test_same_sdf_records_in_a_different_order_keep_the_partition(
    tmp_path, monkeypatch, seed
):
    expected = _assembled(LIBRARY_ORDER_CASE)
    source = _input_smiles(expected[1])
    path = _copy_library(tmp_path, seed)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    result = convert_smiles(source)
    _assert_registered_recovery(expected, result)


def test_registration_changes_unknown_to_known_without_restart(tmp_path, monkeypatch):
    expected_unknown = _assembled(UNKNOWN_CONTEXT)
    source = _input_smiles(expected_unknown[1])
    path = _copy_library(tmp_path)
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    # Warm the converter's caches before adding the previously unseen monomer.
    initial = convert_smiles(source)
    assert Chem.MolToSmiles(_assembled(initial.cabiln)[1]) == Chem.MolToSmiles(
        expected_unknown[1]
    )
    register_monomer("N[C@@H](CCC(F)(F)F)C(=O)O", "AuditNewInternal", sdf_path=path)
    expected_registered = _assembled("A-G-AuditNewInternal-K-V")
    result = convert_smiles(source)
    _assert_registered_recovery(expected_registered, result)


def _assert_assignment_coverage(source, result):
    """Check the public witness in the caller's input index space."""
    atoms = [
        atom for assignment in result.assignments for atom in assignment.source_atoms
    ]
    assert Counter(atoms) == Counter(range(source.GetNumAtoms()))
    for assignment in result.assignments:
        assert isinstance(assignment.symbol, str) and assignment.symbol
        assert isinstance(assignment.recognized, bool)
        assert assignment.source_atoms
        assert len({slot for slot, _ in assignment.attachments}) == len(
            assignment.attachments
        )
        assert all(
            isinstance(slot, int) and slot > 0 and atom in assignment.source_atoms
            for slot, atom in assignment.attachments
        )


@pytest.mark.parametrize("seed", [None, 41])
@pytest.mark.parametrize(
    "smiles,expected_assignments,status",
    [
        (
            "N[C@@H](C)C(=O)NCC(=O)O",
            [
                ("A", range(0, 5), {1: 0, 2: 3, 3: 0}, True),
                ("G", range(5, 10), {1: 5, 2: 7, 3: 5}, True),
            ],
            "complete",
        ),
        (
            "N[C@@H](C)C(=O)NCC(=O)N[C@@H](CCC(F)(F)F)C(=O)"
            "N[C@@H](CCCCN)C(=O)N[C@@H](C(C)C)C(=O)O",
            [
                ("A", range(0, 5), {1: 0, 2: 3, 3: 0}, True),
                ("G", range(5, 9), {1: 5, 2: 7, 3: 5}, True),
                (UNKNOWN, range(9, 19), {1: 9, 2: 17}, False),
                ("K", range(19, 28), {1: 19, 2: 26, 3: 19, 4: 25, 5: 25}, True),
                ("V", range(28, 36), {1: 28, 2: 33, 3: 28}, True),
            ],
            "partial",
        ),
    ],
    ids=["known", "local-unknown"],
)
def test_public_assignments_preserve_original_atoms_and_attachment_sites(
    smiles, expected_assignments, status, seed
):
    # Atom boundaries and anchor indices are authored directly against the
    # explicit source SMILES, independently of recognition and CABILN lowering.
    reference = Chem.MolFromSmiles(smiles)
    value = smiles if seed is None else _input_smiles(reference, seed)
    source = Chem.MolFromSmiles(value)
    result = convert_smiles(value)
    assert result.recognition_status == status
    _assert_assignment_coverage(source, result)
    assert Chem.MolToSmiles(_assembled(result.cabiln)[1]) == Chem.MolToSmiles(source)
    assert len(result.assignments) == len(expected_assignments)
    expected_templates = [
        _template(Sequence(symbol).s_monomers[0])
        for symbol, _, _, _ in expected_assignments
    ]
    actual_templates = [
        _template(Sequence(assignment.symbol).s_monomers[0])
        for assignment in result.assignments
    ]
    assert actual_templates == expected_templates
    mappings = source.GetSubstructMatches(
        reference, useChirality=True, uniquify=False, maxMatches=1024
    )
    assert 0 < len(mappings) < 1024
    assert any(
        all(
            set(assignment.source_atoms) == {mapping[atom] for atom in atoms}
            and dict(assignment.attachments)
            == {slot: mapping[atom] for slot, atom in anchors.items()}
            and assignment.recognized is recognized
            for assignment, (_, atoms, anchors, recognized) in zip(
                result.assignments, expected_assignments
            )
        )
        for mapping in mappings
    ), "Original input atom identities or attachment anchors were changed"


@pytest.mark.parametrize(
    "notation",
    [
        pytest.param("A-G%V", id="disconnected-components"),
        pytest.param("ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am", id="reaction-removes-atoms"),
    ],
)
def test_complete_assignments_use_global_source_indices_without_virtual_atoms(notation):
    expected = _assembled(notation)
    value = _input_smiles(expected[1], 41)
    source = Chem.MolFromSmiles(value)
    result = convert_smiles(value)
    assert result.recognition_status == "complete"
    assert all(assignment.recognized for assignment in result.assignments)
    _assert_assignment_coverage(source, result)
    _assert_registered_recovery(expected, result)
    mappings = source.GetSubstructMatches(
        expected[1], useChirality=True, uniquify=False, maxMatches=1024
    )
    assert 0 < len(mappings) < 1024
    assert any(
        _partition_key([{mapping[atom] for atom in atoms} for atoms in expected[2]])
        == _partition_key(
            [assignment.source_atoms for assignment in result.assignments]
        )
        for mapping in mappings
    ), "Public assignments do not describe the independently constructed residues"


@pytest.mark.parametrize(
    "budgets,warning",
    [
        (RecognitionBudgets(max_states=1), "state limit"),
        (RecognitionBudgets(max_patterns=1), "pattern limit"),
    ],
)
def test_budget_cutoff_is_unresolved_when_only_opaque_cycle_is_available(
    budgets, warning
):
    expected = _assembled("!1-A-G-!1")
    source = _input_smiles(expected[1])
    baseline = convert_smiles(source)
    assert baseline.recognition_status == "complete"
    assert baseline.search_complete
    _assert_registered_recovery(expected, baseline)
    result = convert_smiles(source, budgets=budgets)
    assert result.recognition_status == "unresolved"
    assert not result.search_complete
    assert not result.assignments
    assert not result.details
    assert result.synthetic_components == (0,)
    assert any(warning in message.lower() for message in result.warnings)
    assert any("coarse" in message.lower() for message in result.warnings)
    assert Chem.MolToSmiles(_assembled(result.cabiln)[1]) == Chem.MolToSmiles(
        expected[1]
    )


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param(
            "CC(=O)N[C@H]1Cc2cn(nn2)CCCC[C@@H](C(N)=O)NC(=O)[C@H](C)"
            "NC(=O)[C@H](C)NC1=O",
            id="CuAAC-intramolecular",
        ),
        pytest.param(
            "CC(=O)N[C@@H](CCCCn1cc(C[C@H](NC(C)=O)C(N)=O)nn1)C(N)=O",
            id="CuAAC-intermolecular",
        ),
        pytest.param(
            "CC(=O)N[C@H]1CC2CCCCc3c(nnn3CCCC[C@@H](C(N)=O)NC(=O)[C@H](C)"
            "NC(=O)[C@H](C)NC1=O)C2",
            id="SPAAC-intramolecular",
        ),
        pytest.param(
            "CC(=O)N[C@@]1(C)CCCC=CCCCCCC[C@](C)(C(N)=O)NC(=O)[C@H](C)"
            "NC(=O)[C@H](C)NC1=O",
            id="metathesis",
        ),
    ],
)
def test_reaction_ownership_matches_renderer_residue_map(smiles):
    # Fixed product SMILES keep the molecular/stereo/charge expectation separate
    # from the assembly rule whose provenance is being tested.
    source = Chem.MolFromSmiles(smiles)
    result = convert_smiles(smiles)
    assert result.recognition_status == "complete"
    _assert_assignment_coverage(source, result)
    sequence = Sequence(result.cabiln, warning_sink=lambda _message: None)
    assembly = Molecule(sequence)
    product = assembly.get_molecule("ROMol")
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(source)
    # A nearest-neighbor fallback would hide missing reaction provenance.
    assert all(atom.HasProp("_residue_idx") for atom in product.GetAtoms())
    residue_map = assembly.get_residue_atom_map()
    mappings = source.GetSubstructMatches(
        product, useChirality=True, uniquify=False, maxMatches=1024
    )
    assert 0 < len(mappings) < 1024
    assert any(
        all(
            set(assignment.source_atoms)
            == {mapping[atom] for atom in residue_map[assignment.residue_index]}
            for assignment in result.assignments
        )
        for mapping in mappings
    ), "Renderer ownership differs from the converter's reaction witness"
