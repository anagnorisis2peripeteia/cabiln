"""Shared assembly and atom-partition oracles for chemistry tests."""

from rdkit import Chem

from pyPept.molecule import Molecule
from pyPept.sequence import Sequence

# Sequence uses the bundled CHUCKLES library by default.
_LIB = None


def _assemble(biln):
    seq = Sequence(biln) if _LIB is None else Sequence(biln, monomer_lib=_LIB)
    mol = Molecule(seq)
    return seq, mol


def _romol(biln):
    _, mol = _assemble(biln)
    return mol.get_molecule(fmt='ROMol')


def _smiles(biln):
    return Chem.MolToSmiles(_romol(biln))


def _natoms(biln):
    return _romol(biln).GetNumAtoms()


def _assert_same_monomer_partition(expected, actual):
    """Require the same product and monomer boundaries under graph isomorphism.

    Equivalent aliases, symmetric scaffold arms, and branch serialization can
    differ. Known monomers cannot be replaced by synthetic tokens, even if each
    synthetic token preserves the original monomer's atom boundaries.
    """
    if '<' not in expected:
        assert '<' not in actual, f'Known monomers became synthetic: {actual!r}'
    left_sequence, right_sequence = Sequence(expected), Sequence(actual)
    left = Molecule(left_sequence)
    right = Molecule(right_sequence)
    source, rebuilt = left.get_molecule(fmt='ROMol'), right.get_molecule(fmt='ROMol')
    assert Chem.MolToSmiles(source) == Chem.MolToSmiles(rebuilt)
    source_groups = tuple(map(frozenset, left.get_residue_atom_map().values()))
    rebuilt_groups = set(map(frozenset, right.get_residue_atom_map().values()))
    assert len(source_groups) == len(left_sequence.s_monomers)
    assert len(rebuilt_groups) == len(right_sequence.s_monomers)
    assert set().union(*source_groups) == set(range(source.GetNumAtoms()))
    assert set().union(*rebuilt_groups) == set(range(rebuilt.GetNumAtoms()))
    assert sum(map(len, source_groups)) == source.GetNumAtoms()
    assert sum(map(len, rebuilt_groups)) == rebuilt.GetNumAtoms()
    assert len(source_groups) == len(rebuilt_groups), (expected, actual)
    matches = rebuilt.GetSubstructMatches(
        source, useChirality=True, uniquify=False, maxMatches=4096)
    assert any(
        {frozenset(match[index] for index in group) for group in source_groups}
        == rebuilt_groups for match in matches
    ), f'Monomer atom boundaries changed: {expected!r} -> {actual!r}'
