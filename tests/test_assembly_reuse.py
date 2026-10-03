"""Assembly reuse and batch joins preserve the molecule and its provenance."""

from dataclasses import replace

import pytest
from rdkit import Chem

from pyPept import assembly_cache, bond_plan
from pyPept.interfaces.reaction_library import REACTION_INDEX, _compiled_reaction
from pyPept.molecule import Molecule
from pyPept.peptide import Peptide
from pyPept.sequence import Sequence


@pytest.fixture(autouse=True)
def empty_cache():
    assembly_cache.clear_assembly_cache()
    yield
    assembly_cache.clear_assembly_cache()


def lineage(assembly):
    molecule = Chem.Mol(assembly.mol)
    for atom in molecule.GetAtoms():
        original = atom.GetIntProp('_template_atom') if atom.HasProp('_template_atom') else 999
        atom.SetAtomMapNum(1 + 1000 * atom.GetIntProp('_residue_idx') + original)
    bonds = sorted(
        (bond.GetIntProp('_connection_idx'), tuple(sorted((
            bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        ))))
        for bond in molecule.GetBonds() if bond.HasProp('_connection_idx')
    )
    return Chem.MolToSmiles(molecule), bonds


@pytest.mark.parametrize('source', [
    'ac-A-DVal-P-G-am',
    'C.!r(4,4)-C.!r',
    'K.!r(2,4).!r(4,2)',
    'K.[G(4,2).ac(1,2)]-A',
    'K.ac(4,2).ac(5,2)',
    'K._Me(4,2)._Me(5,2)',
    'G-AzK-G',
    'ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TBMB.!1.!2.!3',
])
def test_batch_join_preserves_general_executor_product_stereo_and_lineage(source, monkeypatch):
    batch = bond_plan.batch_join
    joined = []

    def record(*args):
        result = batch(*args)
        joined.append(result is not None)
        return result

    monkeypatch.setattr(bond_plan, 'batch_join', record)
    fast = Molecule(Sequence(source), depiction=None)
    assert joined == [True]
    assembly_cache.clear_assembly_cache()
    monkeypatch.setattr(bond_plan, 'batch_join', lambda *args: None)
    general = Molecule(Sequence(source), depiction=None)
    assert Chem.MolToSmiles(fast.mol) == Chem.MolToSmiles(general.mol)
    assert lineage(fast) == lineage(general)
    for endpoint, index in fast.get_attachment_atom_map().items():
        other = general.mol.GetAtomWithIdx(general.get_attachment_atom_map()[endpoint])
        atom = fast.mol.GetAtomWithIdx(index)
        assert (atom.GetAtomicNum(), atom.GetIsotope()) == (other.GetAtomicNum(), other.GetIsotope())


@pytest.mark.parametrize('source', [
    'Pra.!r(4,4)-A-AzK.!r',
    'G.<[1*]/C=C/C>_(3,1)',
    'A-<[1*]N[C@@H](C[O-])C([2*])=O.[Na+]>-G',
])
def test_complex_reactions_and_nonlocal_port_edits_use_general_executor(source, monkeypatch):
    batch = bond_plan.batch_join
    outcomes = []

    def record(*args):
        result = batch(*args)
        outcomes.append(result is None)
        return result

    monkeypatch.setattr(bond_plan, 'batch_join', record)
    assembly = Molecule(Sequence(source), depiction=None)
    assert outcomes == [True]
    assembly_cache.clear_assembly_cache()
    monkeypatch.setattr(bond_plan, 'batch_join', lambda *args: None)
    general = Molecule(Sequence(source), depiction=None)
    assert lineage(assembly) == lineage(general)


def test_cached_join_queries_survive_reaction_cache_eviction():
    smirks = REACTION_INDEX['backbone_n', 'backbone_c']['steps'][0]
    nitrogen, carbon = bond_plan._join_queries(smirks)
    _compiled_reaction.cache_clear()
    assert bond_plan._matches(Chem.MolFromSmiles('C[5*]'), nitrogen, 5) is False
    assert bond_plan._matches(Chem.MolFromSmiles('N[5*]'), nitrogen, 5)
    assert bond_plan._matches(Chem.MolFromSmiles('CC(=O)[6*]'), carbon, 6)


def test_reuse_detaches_drawings_and_retains_ports_origins_and_connections(monkeypatch):
    source = 'A-G'
    first = Molecule(Sequence(source), depiction=None)
    expected = lineage(first)
    sites = first.current_sites(1)
    first.mol.GetAtomWithIdx(0).SetIsotope(15)
    first._port_mol.GetAtomWithIdx(0).SetIntProp('_template_atom', 900)
    first._port_labels.clear()

    def forbidden(*args):
        raise AssertionError('An unchanged resolved graph should reuse assembly')

    monkeypatch.setattr(Molecule, '_Molecule__assemble', forbidden)
    reused = Molecule(Sequence(source), depiction='rdkit')
    assert lineage(reused) == expected
    assert reused.current_sites(1) == sites
    assert reused.mol.GetNumConformers() == 1


def test_reuse_is_bound_to_template_content_and_rules(monkeypatch):
    peptide = Peptide.from_sequence(Sequence('G'))
    node = peptide.occurrences[0]
    first = Molecule(peptide, depiction=None)
    template = node.definition.copy_template()
    template.GetAtomWithIdx(1).SetIsotope(15)
    edited = replace(peptide, occurrences=(replace(node, definition=replace(node.definition, _template=template)),))
    second = Molecule(edited, depiction=None)
    assert Chem.MolToSmiles(first.mol) != Chem.MolToSmiles(second.mol)
    original_key = assembly_cache.assembly_key(peptide)
    monkeypatch.setattr(assembly_cache, 'REACTION_FINGERPRINT', b'changed rules')
    assert assembly_cache.assembly_key(peptide) != original_key


def test_reuse_evicts_by_bytes_and_recency(monkeypatch):
    def build(source):
        return Molecule(Sequence(source), depiction=None)

    build('G')
    one_size = assembly_cache._bytes
    monkeypatch.setattr(assembly_cache, '_MAX_BYTES', one_size * 2)
    build('A')
    build('G')
    key = assembly_cache.assembly_key(Peptide.from_sequence(Sequence('G')))
    assert assembly_cache.get_assembly(key) is not None
    build('A-G-G-A')
    assert assembly_cache._bytes <= one_size * 2
    monkeypatch.setattr(assembly_cache, '_MAX_BYTES', 32 * 1024 * 1024)
    monkeypatch.setattr(assembly_cache, '_MAX_ENTRIES', 2)
    build('P')
    build('V')
    assert assembly_cache.get_assembly(key) is None
