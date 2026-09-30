"""Recognition must expose useful known pieces, not just lossless encoding."""

from collections import Counter
from dataclasses import replace

import pytest
from rdkit import Chem

from pyPept import monomer_store
from pyPept.interfaces.cli_monomer import register_monomer
from pyPept.molecule import Molecule
from pyPept.recognition import (
    RecognitionBudgets,
    RejectedRegion,
    _compile_patterns,
    _enumerate_candidates,
    _patterns,
    _template_core,
    recognize,
)
from pyPept.recognition_reactions import recognition_variants
from pyPept.sequence import Sequence


@pytest.fixture
def select_library(tmp_path, monkeypatch):
    _, bundled = monomer_store._load_sdf()

    def select(*symbols):
        path = tmp_path / "selected.sdf"
        with Chem.SDWriter(str(path)) as writer:
            for symbol in symbols:
                writer.write(bundled[symbol])
        monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
        return path

    return select


def _recognize(smiles, **budget_overrides):
    return recognize(
        Chem.MolFromSmiles(smiles), budgets=RecognitionBudgets(**budget_overrides)
    )


def _best(result):
    assert result.covers
    return result.covers[0]


def _assert_partition(molecule, cover, *, complete):
    owned = set()
    owners = {}
    for candidate in cover:
        assert not owned & candidate.atoms
        owned.update(candidate.atoms)
        for atom in candidate.atoms:
            owners[atom] = candidate
        for inside, outside, slot in candidate.ports:
            assert inside in candidate.atoms
            assert outside not in candidate.atoms
            assert molecule.GetBondBetweenAtoms(inside, outside) is not None
            assert (inside, slot) in candidate.anchors
    for candidate in cover:
        for inside, outside, _ in candidate.ports:
            if outside in owners:
                assert any(
                    (a, b) == (outside, inside) for a, b, _ in owners[outside].ports
                )
    if complete:
        assert owned == set(range(molecule.GetNumAtoms()))
    return owned


def test_known_dipeptide_has_exact_ownership_and_all_attachment_anchors(select_library):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
    result = recognize(source)
    cover = _best(result)
    _assert_partition(source, cover, complete=True)
    by_name = {candidate.symbol: candidate for candidate in cover}
    assert by_name["G"].atoms == frozenset((0, 1, 2, 3))
    assert by_name["A"].atoms == frozenset(range(4, 10))
    assert by_name["G"].ports == ((2, 4, 2),)
    assert by_name["A"].ports == ((4, 2, 1),)
    assert by_name["G"].anchors == ((0, 1), (0, 3), (2, 2))
    assert dict(by_name["G"].attachment_types)[3] == "backbone_n_mod"
    assert not result.exhausted and not result.match_truncated


def test_unknown_middle_keeps_both_known_neighbors_and_real_boundary_ports(
    select_library,
):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](CCC(F)(F)F)C(=O)N[C@@H](C)C(=O)O")
    cover = _best(recognize(source))
    owned = _assert_partition(source, cover, complete=False)
    by_name = {candidate.symbol: candidate for candidate in cover}
    assert set(by_name) == {"G", "A"}
    assert by_name["G"].atoms == frozenset(range(4))
    assert by_name["A"].atoms == frozenset(range(14, 20))
    assert by_name["G"].ports == ((2, 4, 2),)
    assert by_name["A"].ports == ((14, 12, 1),)
    assert set(range(source.GetNumAtoms())) - owned == set(range(4, 14))


@pytest.mark.parametrize(
    "smiles",
    ["[NH3+]CC(=O)O", "N[13CH2]C(=O)O"],
    ids=["charged-amine", "isotope"],
)
def test_neutral_unlabeled_template_does_not_claim_changed_atoms(
    select_library, smiles
):
    select_library("G")
    assert _best(_recognize(smiles)) == ()


def test_charged_leaving_group_is_unknown_instead_of_claimed_as_neutral(select_library):
    select_library("G")
    cover = _best(_recognize("NCC(=O)[O-]"))
    assert len(cover) == 1 and cover[0].symbol == "G"
    # The glycine core is exact, but its terminal O- cannot be owned as the
    # library's restored OH. A caller must preserve and validate that boundary.
    assert cover[0].atoms == frozenset(range(4))
    assert cover[0].ports == ((2, 4, 2),)


def test_defined_stereo_is_preserved_and_unspecified_inference_is_explicit_policy(
    select_library,
):
    select_library("A")
    assert [item.symbol for item in _best(_recognize("N[C@@H](C)C(=O)O"))] == ["A"]
    assert _best(_recognize("N[C@H](C)C(=O)O")) == ()
    assert [item.symbol for item in _best(_recognize("NC(C)C(=O)O"))] == ["A"]
    assert _best(_recognize("NC(C)C(=O)O", allow_unspecified_stereo=False)) == ()


def test_unspecified_template_cannot_erase_defined_source_stereo(select_library):
    path = select_library("A")
    molecule = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
    Chem.RemoveStereochemistry(molecule)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule)
    assert _best(_recognize("N[C@@H](C)C(=O)O")) == ()


def test_defined_double_bond_stereo_cannot_match_opposite_or_unspecified_template(
    select_library,
):
    path = select_library("G")
    trans = "N[C@@H](C/C=C/F)C(=O)O"
    cis = "N[C@@H](C/C=C\\F)C(=O)O"
    register_monomer(trans, "FluoroAlkene")
    assert {item.symbol for item in _best(_recognize(trans))} == {"FluoroAlkene"}
    assert "FluoroAlkene" not in {item.symbol for item in _best(_recognize(cis))}
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    for bond in molecules[-1].GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    with Chem.SDWriter(str(path)) as writer:
        for molecule in molecules:
            writer.write(molecule)
    assert "FluoroAlkene" not in {item.symbol for item in _best(_recognize(trans))}


def test_distinct_nitrogen_slots_are_not_deduplicated(select_library):
    select_library("G", "A")
    result = _recognize("NCC(=O)N[C@@H](C)C(=O)O")
    alanine_slots = {
        next(item for item in cover if item.symbol == "A").ports[0][2]
        for cover in result.covers
    }
    assert alanine_slots == {1, 3}
    assert next(item for item in _best(result) if item.symbol == "A").ports[0][2] == 1


def test_backbone_exit_alternatives_are_ranked_by_supported_r2_r1_connections(
    select_library,
):
    select_library("ac", "am", "G", "E", "E_g")
    source = Molecule(Sequence("ac-E_g-G-am")).get_molecule(fmt="ROMol")
    cover = _best(recognize(source))
    _assert_partition(source, cover, complete=True)
    assert {candidate.symbol for candidate in cover} == {"ac", "am", "G", "E_g"}


def test_equivalent_name_choice_is_deterministic_without_a_name_priority_table(
    select_library,
):
    path = select_library("G")
    original = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
    alias = Chem.Mol(original)
    alias.SetProp("symbol", "LongGlycineName")
    alias.SetProp("m_abbr", "LongGlycineName")
    with Chem.SDWriter(str(path)) as writer:
        writer.write(alias)
        writer.write(original)
    result = _recognize("NCC(=O)O")
    assert _best(result)[0].symbol == "G"
    assert {cover[0].symbol for cover in result.covers} == {"G", "LongGlycineName"}


def test_missing_declared_types_uses_shared_effective_detection(select_library):
    path = select_library("G")
    molecule = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
    molecule.ClearProp("m_chem_types")
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule)
    candidate = _best(_recognize("NCC(=O)O"))[0]
    assert candidate.symbol == "G" and candidate.has_backbone
    assert dict(candidate.attachment_types)[1] == "backbone_n"


def test_library_registration_changes_recognition_without_restart(select_library):
    select_library("G", "A")
    source = "NCC(=O)N[C@@H](CCC(F)(F)F)C(=O)N[C@@H](C)C(=O)O"
    assert {item.symbol for item in _best(_recognize(source))} == {"G", "A"}
    register_monomer("N[C@@H](CCC(F)(F)F)C(=O)O", "NewFluoro")
    molecule = Chem.MolFromSmiles(source)
    cover = _best(recognize(molecule))
    assert {item.symbol for item in cover} == {"G", "A", "NewFluoro"}
    _assert_partition(molecule, cover, complete=True)


@pytest.mark.parametrize(
    "notation,symbols,expected",
    [
        (
            "ac-A-redG2-G-am",
            ("ac", "A", "redG2", "G", "am"),
            {"ac", "A", "redG2", "G", "am"},
        ),
        ("!1-A-G-C-A-!1", ("A", "G", "C"), {"A", "G", "C"}),
        (
            "ac-K.[G(4,2).ac(1,2)]-C.!1(4,4)-am%ac-C.!1-G-am",
            ("ac", "K", "G", "C", "am"),
            {"ac", "K", "G", "C", "am"},
        ),
    ],
    ids=["mixed-backbone", "cycle", "branched-two-chain-disulfide"],
)
def test_supported_topologies_have_complete_known_partitions(
    select_library, notation, symbols, expected
):
    select_library(*symbols)
    source = Molecule(Sequence(notation, warning_sink=lambda _: None)).get_molecule(
        fmt="ROMol"
    )
    result = recognize(source)
    cover = _best(result)
    assert {candidate.symbol for candidate in cover} == expected
    _assert_partition(source, cover, complete=True)


def test_state_limit_returns_observed_partial_progress_and_explicit_status(
    select_library,
):
    select_library("G", "A")
    result = _recognize("NCC(=O)N[C@@H](C)C(=O)O", max_states=2)
    assert result.exhausted and result.states == 2
    assert len(_best(result)) == 1
    assert any("Search state limit" in warning for warning in result.warnings)


def test_match_limit_is_explicit_not_a_claim_of_exhaustive_recognition(select_library):
    select_library("G")
    result = _recognize("NCC(=O)O.NCC(=O)O", max_matches_per_pattern=1)
    assert result.match_truncated
    assert any("Match limit" in warning for warning in result.warnings)


def test_alternative_limit_keeps_preferred_backbone_slots_but_reports_ambiguity(
    select_library,
):
    select_library("G", "A")
    result = _recognize("NCC(=O)N[C@@H](C)C(=O)O", max_covers=1)
    assert result.exhausted
    assert next(item for item in _best(result) if item.symbol == "A").ports[0][2] == 1
    assert any("Alternative limit" in warning for warning in result.warnings)


def test_examined_solution_limit_is_explicit_and_keeps_observed_preferred_slots(
    select_library,
):
    select_library("G", "A")
    result = _recognize("NCC(=O)N[C@@H](C)C(=O)O", max_solutions=1)
    assert result.exhausted
    assert next(item for item in _best(result) if item.symbol == "A").ports[0][2] == 1
    assert any("Examined solution limit" in warning for warning in result.warnings)


def test_solution_limit_does_not_starve_a_better_known_atom_partition(select_library):
    path = select_library("G")
    with Chem.SDWriter(str(path)) as writer:
        for name, chain in (("Four", "CCCC"), ("Three", "CCC"), ("Two", "CC")):
            template = Chem.MolFromSmiles(f"[1*]{chain}[2*]")
            template.SetProp("symbol", name)
            template.SetProp("m_abbr", name)
            template.SetProp("m_type", "linker")
            template.SetProp("m_Rgroups", "[H],[H]")
            writer.write(template)
    source = Chem.MolFromSmiles("CCCCC")
    result = recognize(source, budgets=RecognitionBudgets(max_solutions=1))
    # A tempting four-carbon match leaves one carbon unknown. Its first partial
    # solution must not consume the limit before the 3+2 complete partition.
    cover = _best(result)
    _assert_partition(source, cover, complete=True)
    assert sorted(len(candidate.atoms) for candidate in cover) == [2, 3]


def test_rejected_high_coverage_covers_do_not_prune_admissible_partial_covers(
    select_library,
):
    select_library("G", "A")
    calls = Counter()

    def admissible(cover):
        calls[cover] += 1
        # Model a region whose apparent library match cannot be emitted without
        # changing chemistry. The other known residue must remain available.
        return all(candidate.symbol != "A" for candidate in cover)

    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
    result = recognize(
        source,
        budgets=RecognitionBudgets(max_solutions=1),
        accept_cover=admissible,
    )
    cover = _best(result)
    assert len(cover) == 1 and cover[0].symbol == "G"
    assert cover[0].atoms == frozenset(range(4))
    assert any(
        any(candidate.symbol == "A" for candidate in checked) for checked in calls
    )
    assert all(count == 1 for count in calls.values())


def test_admissibility_programming_errors_are_not_swallowed(select_library):
    select_library("G")

    def invalid_callback(_cover):
        raise TypeError("broken validation callback")

    with pytest.raises(TypeError, match="broken validation callback"):
        recognize(Chem.MolFromSmiles("NCC(=O)O"), accept_cover=invalid_callback)


def test_alias_products_do_not_starve_the_required_unknown_boundary(select_library):
    path = select_library("am", "G", "A")
    entries = list(Chem.SDMolSupplier(str(path), removeHs=False))
    alias = Chem.Mol(next(mol for mol in entries if mol.GetProp("symbol") == "G"))
    alias.SetProp("symbol", "AlternateGlycine")
    alias.SetProp("m_abbr", "AlternateGlycine")
    with Chem.SDWriter(str(path)) as writer:
        for entry in entries + [alias]:
            writer.write(entry)
    # A source-specific histidine tautomer precedes twenty independently written
    # known residues. Its free alpha N superficially matches an amine cap.
    source = Chem.MolFromSmiles(
        "N[C@@H](Cc1c[nH]cn1)C(=O)" + "NCC(=O)N[C@@H](C)C(=O)" * 10 + "O"
    )
    calls = []

    def boundary_admissible(cover):
        calls.append(cover)
        # The emitter requires this N to remain with the unknown residue. This
        # models the local boundary constraint without filtering a monomer name.
        return not any(0 in candidate.atoms for candidate in cover)

    result = recognize(
        source,
        budgets=RecognitionBudgets(max_states=500, max_assignments_per_partition=4),
        accept_cover=boundary_admissible,
    )
    cover = _best(result)
    assert len(cover) == 20
    assert all(candidate.has_backbone for candidate in cover)
    assert not any(0 in candidate.atoms for candidate in cover)
    assert any(any(0 in candidate.atoms for candidate in checked) for checked in calls)
    assert len(calls) <= 16
    assert result.states < 500
    assert result.exhausted
    assert any("Slot/name assignment limit" in warning for warning in result.warnings)


def test_grouped_slot_assignments_retain_alternative_backbone_exits(select_library):
    select_library("ac", "am", "G", "E", "E_g")
    source = Molecule(Sequence("ac-E_g-G-am")).get_molecule(fmt="ROMol")
    result = recognize(
        source,
        budgets=RecognitionBudgets(max_assignments_per_partition=128),
        accept_cover=lambda _cover: True,
    )
    _assert_partition(source, _best(result), complete=True)
    assert {candidate.symbol for candidate in _best(result)} == {"ac", "am", "G", "E_g"}
    assert any(
        any(candidate.symbol == "E" for candidate in cover) for cover in result.covers
    )


def test_grouped_assignment_budget_does_not_claim_exhaustive_slot_validation(
    select_library,
):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
    result = recognize(
        source,
        budgets=RecognitionBudgets(max_assignments_per_partition=1),
        accept_cover=lambda _cover: True,
    )
    _assert_partition(source, _best(result), complete=True)
    assert result.exhausted
    assert any(
        "Slot/name assignment limit (1" in warning for warning in result.warnings
    )


def _mixed_connection_problem(*, max_states=100):
    from pyPept.attachments import reaction_for_types
    from pyPept.recognition import Candidate, _RecognitionProblem

    # Eight early names have the same ownership but incompatible site chemistry.
    # Both later alkyl-halide aliases must remain eligible for the thiol partner.
    choices = tuple(
        Candidate(
            name, frozenset({0}), ((0, 1, 4),), ((0, 4),), False, "chem", 1,
            ((4, chemistry),),
        )
        for name, chemistry in (
            *((f"A{index}", "carbon") for index in range(8)),
            ("ZValid0", "alkyl_halide_c"), ("ZValid1", "alkyl_halide_c"),
        )
    )
    partner = Candidate(
        "Cap", frozenset({1}), ((1, 0, 4),), ((1, 4),), False, "chem", 1,
        ((4, "thiol"),),
    )
    return _RecognitionProblem(
        Chem.MolFromSmiles("CS"),
        RecognitionBudgets(max_states=max_states, max_assignments_per_partition=2),
        (*choices, partner), False, (), connection_policy=reaction_for_types,
    )


def test_connection_policy_keeps_later_supported_aliases_in_ownership_group():
    problem = _mixed_connection_problem()
    result = problem.search(lambda cover: sum(len(c.atoms) for c in cover) == 2)
    assert {c.symbol for cover in result.covers for c in cover} == {
        "ZValid0", "ZValid1", "Cap",
    }
    for cover in result.covers:
        _assert_partition(problem.molecule, cover, complete=True)
    # Raw recognition retains its structural proposal contract when no policy
    # is requested, even for connections outside the reaction library.
    structural = replace(problem, connection_policy=None).search()
    assert any(c.symbol.startswith("A") for cover in structural.covers for c in cover)


def test_connection_policy_rejections_consume_existing_state_budget():
    problem = _mixed_connection_problem(max_states=8)
    result = problem.search(lambda cover: sum(len(c.atoms) for c in cover) == 2)
    assert result.states == 8
    assert not result.covers
    assert result.exhausted
    assert any("Search state limit (8)" in message for message in result.warnings)


def test_hereditary_unknown_region_proof_prunes_a_partial_sidechain_match(
    select_library,
):
    select_library("G", "A", "Orn", "am")
    source = Chem.MolFromSmiles(
        "NCC(=O)N[C@@H](CCCN/C(N)=N/[H])C(=O)" + "NCC(=O)N[C@@H](C)C(=O)" * 6 + "O"
    )
    rejections = []

    def admissible(cover):
        unknown = set(range(source.GetNumAtoms())) - {
            atom for candidate in cover for atom in candidate.atoms
        }
        while unknown:
            region, frontier = set(), [min(unknown)]
            while frontier:
                atom = frontier.pop()
                if atom not in unknown:
                    continue
                unknown.remove(atom)
                region.add(atom)
                frontier.extend(
                    neighbor.GetIdx()
                    for neighbor in source.GetAtomWithIdx(atom).GetNeighbors()
                    if neighbor.GetIdx() in unknown
                )
            has_carbonyl = any(
                source.GetAtomWithIdx(atom).GetAtomicNum() == 6
                and any(
                    bond.GetBondType() == Chem.BondType.DOUBLE
                    and bond.GetOtherAtomIdx(atom) in region
                    and bond.GetOtherAtom(source.GetAtomWithIdx(atom)).GetAtomicNum()
                    == 8
                    for bond in source.GetAtomWithIdx(atom).GetBonds()
                )
                for atom in region
            )
            if not has_carbonyl:
                rejections.append(region)
                return RejectedRegion(
                    frozenset(region),
                    "Every unknown piece requires a carbonyl",
                    requires_complete_recognition=True,
                )
        return True

    result = recognize(
        source,
        budgets=RecognitionBudgets(max_states=1000),
        accept_cover=admissible,
    )
    cover = _best(result)
    assert all(candidate.symbol != "Orn" for candidate in cover)
    assert sum(candidate.symbol in ("G", "A") for candidate in cover) == 13
    assert 0 < len(rejections) < 20
    assert result.states < 1000


def test_region_with_a_possible_recognized_piece_does_not_create_a_conflict(
    select_library,
):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")

    def admissible(cover):
        if len(cover) != 1 or cover[0].symbol != "A":
            return False
        if cover[0].ports[0][2] == 1:
            # Glycine could fill this entire region while A stays selected.
            # The engine must decline this insufficient conflict certificate.
            return RejectedRegion(frozenset(range(4)), "Rejected proposal")
        return True

    cover = _best(recognize(source, accept_cover=admissible))
    assert len(cover) == 1 and cover[0].symbol == "A"
    assert cover[0].ports[0][2] == 3


def test_rejection_cannot_learn_from_atoms_already_owned_by_the_cover(select_library):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")

    def admissible(cover):
        if len(cover) == 2:
            return RejectedRegion(frozenset((0,)), "Invalid owned-region evidence")
        return len(cover) == 1 and cover[0].symbol == "G"

    cover = _best(recognize(source, accept_cover=admissible))
    assert len(cover) == 1 and cover[0].symbol == "G"


def test_frontier_limit_reports_discarded_work_explicitly(select_library):
    select_library("G", "A")
    source = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
    result = recognize(
        source,
        budgets=RecognitionBudgets(max_frontier_states=1),
        accept_cover=lambda _cover: True,
    )
    _assert_partition(source, _best(result), complete=True)
    assert result.exhausted
    assert any("Ownership frontier limit (1)" in warning for warning in result.warnings)


def test_repeated_twenty_residue_chain_does_not_enumerate_every_nitrogen_slot_assignment(
    select_library,
):
    select_library("G", "A")
    source = Molecule(Sequence("-".join(["G", "A"] * 10))).get_molecule(fmt="ROMol")
    result = recognize(source)
    cover = _best(result)
    _assert_partition(source, cover, complete=True)
    assert len(cover) == 20
    owners = {atom: candidate for candidate in cover for atom in candidate.atoms}
    backbone_edges = sum(
        slot == 2 and (outside, inside, 1) in owners[outside].ports
        for candidate in cover
        for inside, outside, slot in candidate.ports
    )
    assert backbone_edges == 19
    assert result.states < 5000
    assert result.exhausted  # The omitted alternatives remain disclosed.


def test_pattern_limit_reports_incomplete_template_states(select_library):
    select_library("G")
    result = _recognize("NCC(=O)O", max_patterns=1)
    assert result.match_truncated
    assert any("Pattern limit" in warning for warning in result.warnings)


def test_irrelevant_template_cores_do_not_consume_the_state_expansion_budget(
    select_library,
):
    select_library("G", "A", "W", "K")
    source = Chem.MolFromSmiles("NCC(=O)O")
    # Glycine has eight attachment masks. None of the larger cores can occur in
    # this input, so their masks must not consume compilation time or memory.
    result = recognize(source, budgets=RecognitionBudgets(max_patterns=8))
    assert {candidate.symbol for candidate in _best(result)} == {"G"}
    _assert_partition(source, _best(result), complete=True)
    assert not result.match_truncated


@pytest.mark.parametrize(
    "template,source",
    [
        ("[1*][N+](C)([2H])[13C@@H](F)C([2*])=O", "CN[C@H](F)CO"),
        ("[1*]c1ccccc1", "C1CCCCC1"),
        ("[1*]C([2*])([3*])C([4*])=O", "CCO"),
    ],
    ids=["charge-isotope-stereo-explicit-H", "aromatic-bond-state", "slot-degree"],
)
def test_necessary_core_ignores_state_specific_chemistry(template, source):
    query = _template_core(Chem.MolFromSmiles(template))
    assert Chem.MolFromSmiles(source).HasSubstructMatch(query, useChirality=False)


def test_necessary_core_retains_required_elements_and_ring_topology():
    query = _template_core(Chem.MolFromSmiles("[1*]c1ccncc1"))
    assert not Chem.MolFromSmiles("C1CCCCC1").HasSubstructMatch(query)
    assert not Chem.MolFromSmiles("CCCNCC").HasSubstructMatch(query)
    assert Chem.MolFromSmiles("C1CCNCC1").HasSubstructMatch(query)


@pytest.mark.parametrize(
    "notation",
    [
        "ac-A-redG2-G-am",
        "!1-A-G-C-A-!1",
        "ac-K.[G(4,2).ac(1,2)]-C.!1(4,4)-am%ac-C.!1-G-am",
        "G-<N[C@@H](CCC(F)(F)F)C(=O)O>-A",
        "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am",
        "ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am",
    ],
    ids=["mixed-backbone", "cycle", "branch-disulfide", "unknown", "CuAAC", "RCM"],
)
def test_core_prefilter_preserves_every_exact_candidate_in_all_proposals(
    select_library, notation
):
    select_library(
        "ac", "am", "G", "A", "redG2", "K", "C", "Pra", "AzK", "S5", "R8", "W"
    )
    source = Molecule(Sequence(notation, warning_sink=lambda _: None)).get_molecule(
        fmt="ROMol"
    )
    molecules, _ = monomer_store._load_sdf()
    budgets = RecognitionBudgets()
    unfiltered, truncated, _ = _compile_patterns(molecules, budgets.max_patterns)
    assert not truncated
    for variant in recognition_variants(source).variants:
        filtered, truncated, _ = _patterns(budgets, variant.molecule)
        assert not truncated
        actual = _enumerate_candidates(variant.molecule, filtered, budgets)
        expected = _enumerate_candidates(variant.molecule, unfiltered, budgets)
        assert not actual[1] and not expected[1]
        assert actual[0] == expected[0]


def test_invalid_budgets_fail_clearly():
    with pytest.raises(ValueError, match="max_states"):
        replace(RecognitionBudgets(), max_states=0)


def test_candidate_and_search_do_not_mutate_source(select_library):
    select_library("G", "A")
    molecule = Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)O")
    original = Chem.MolToMolBlock(molecule)
    first = recognize(molecule)
    second = recognize(molecule)
    assert first == second
    assert Chem.MolToMolBlock(molecule) == original
