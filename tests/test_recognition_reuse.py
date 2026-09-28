"""Expensive recognition facts are local to one unchanged conversion input."""

from collections import Counter

from rdkit import Chem

from pyPept import recognition, smiles
from pyPept.structure import compare_structures


def test_partial_import_reuses_candidates_and_verification_across_search_phases(
    monkeypatch,
):
    source = "NCC(=O)N[C@@H](CO[13CH3])C(=O)NCC(=O)O"
    enumeration, assemblies, phases = Counter(), Counter(), []
    enumerate_candidates = recognition._enumerate_candidates
    assemble = smiles._assemble
    search = recognition._RecognitionProblem.search

    def counted_candidates(molecule, *args):
        enumeration[id(molecule)] += 1
        return enumerate_candidates(molecule, *args)

    def counted_assembly(notation, *args):
        assemblies[notation] += 1
        return assemble(notation, *args)

    def counted_search(problem, accept_cover=None):
        phases.append(accept_cover is not None)
        return search(problem, accept_cover)

    monkeypatch.setattr(recognition, "_enumerate_candidates", counted_candidates)
    monkeypatch.setattr(smiles, "_assemble", counted_assembly)
    monkeypatch.setattr(recognition._RecognitionProblem, "search", counted_search)
    result = smiles.convert_smiles(source)
    assert result.recognition_status == "partial"
    assert False in phases and True in phases
    assert all(count == 1 for count in enumeration.values())
    assert all(count == 1 for count in assemblies.values())
    assert len(result.assignments) == 3
    assert compare_structures(Chem.MolFromSmiles(source), assemble(result.cabiln)).exact
