"""One occurrence/site model must survive layout changes and real assembly."""

from dataclasses import replace

import pytest
from rdkit import Chem

from pyPept.molecule import Molecule
from pyPept.peptide import AtomProvenance, Connection, Endpoint, Peptide, serialize
from pyPept.sequence import Sequence
from pyPept.smiles import convert_smiles


def _sequence(source):
    return Sequence(source, track_source=True, warning_sink=lambda _: None)


def _molecule(sequence):
    return Molecule(sequence, depiction=None).get_molecule(fmt="ROMol")


def _assert_serialization(sequence, peptide, notation):
    emitted = serialize(peptide, notation=notation)
    parsed = _sequence(emitted.text)
    restored = Peptide.from_sequence(parsed, emitted.occurrence_order)
    assert set(restored.connections) == set(peptide.connections)
    assert {item.id: item.symbol for item in restored.occurrences} == {
        item.id: item.symbol for item in peptide.occurrences
    }
    assert Chem.MolToSmiles(_molecule(sequence)) == Chem.MolToSmiles(_molecule(parsed))
    return emitted


SOURCES = [
    "ac-A-G-am",
    "A%K.{G(4,2)}-A",
    "ac-K.[G(4,2)]-K.{A(4,2)}-am",
    "ac-K.{G(4,2)[.A(1,2)]}-D-am",
    "ac-K.{G(4,2)[.A(1,2).!1(1,4)]}-D.!1(4,1)-am",
    "ac-K.[G(4,2)[.!1(1,4)]]-D.!1(4,1)-am",
    "ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am",
    "ac-K.G(4,2)-am%G",
    "!1-C.!2(4,4)-A-C.!2-G-!1",
    "ac-C.[TBMB(4,4).!2(5,4).!3(6,4)]-C.!2-C.!3-am",
    "ac-K.[G(4,2).A(1,2)]-am",
    "ac-K.[K(4,2)[.G(1,2)].!1(4,2)]-am%G.!1(2,4)",
    "A-<[1*]N[C@@H](C[O-])C([2*])=O>-G",
    "A-<[1*]N[C@@H](CC%10CCCCC%10)C([2*])=O>-G",
    "  Ala-Gly  %  K-A  ",
]


@pytest.mark.parametrize("source", SOURCES)
def test_preserve_returns_original_bytes_and_occurrence_order(source):
    sequence = _sequence(source)
    peptide = Peptide.from_sequence(sequence)
    emission = serialize(peptide, notation="preserve")
    assert emission.text.encode() == source.encode()
    assert emission.occurrence_order == tuple(range(len(sequence.s_monomers)))


@pytest.mark.parametrize("notation", ["percent", "bracket"])
@pytest.mark.parametrize("source", SOURCES)
def test_layout_changes_preserve_every_occurrence_slot_and_assembled_structure(
    source, notation
):
    sequence = _sequence(source)
    peptide = Peptide.from_sequence(sequence)
    emitted = _assert_serialization(sequence, peptide, notation)
    if ".{" in source:
        for group in peptide.layout.groups:
            if group.protected:
                assert source[group.span.start : group.span.end] in emitted.text


def test_layout_records_branches_on_every_explicit_segment():
    peptide = Peptide.from_sequence(_sequence("A%K.{G(4,2)}-A"))
    assert [segment.roots for segment in peptide.layout.segments] == [(0,), (1, 2)]
    assert peptide.layout.segments[1].members == (1, 2, 3)
    group = peptide.layout.groups[0]
    assert group.host == 1 and group.members == (3,) and group.protected
    assert (group.opening, group.closing) == ("{", "}")


def test_layout_keeps_mixed_delimiters_and_nested_arm_host_identity():
    mixed = Peptide.from_sequence(_sequence("ac-K.[G(4,2)]-K.{A(4,2)}-am"))
    assert [
        (group.host, group.opening, group.protected) for group in mixed.layout.groups
    ] == [(1, "[", False), (2, "{", True)]
    nested = Peptide.from_sequence(_sequence("ac-K.{G(4,2)[.A(1,2)]}-D-am"))
    outer, arm = nested.layout.groups
    assert (outer.host, outer.members, outer.parent) == (1, (4, 5), None)
    assert (arm.host, arm.members, arm.parent) == (4, (5,), outer.id)


def test_legacy_nested_groups_publish_real_semantic_delimiters():
    source = "ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am"
    peptide = Peptide.from_sequence(_sequence(source))
    assert [(group.opening, group.closing) for group in peptide.layout.groups] == [
        ("[", "]"),
        ("[", "]"),
        ("[", "]"),
    ]
    assert serialize(peptide, "preserve").text == source


@pytest.mark.parametrize("prefix", ["", "A%"])
def test_nested_marker_belongs_to_the_arm_that_displays_its_owner(prefix):
    source = prefix + "ac-K.{G(4,2)[.A(1,2).!1(1,4)]}-D.!1(4,1)-am"
    peptide = Peptide.from_sequence(_sequence(source))
    outer, arm = peptide.layout.groups
    nested, explicit = peptide.layout.markers
    assert arm.parent == outer.id
    assert nested.endpoint.occurrence_id in arm.members
    assert nested.group == arm.id
    assert source[nested.span.start : nested.span.end] == ".!1(1,4)"
    assert explicit.group is None


def test_marker_only_nested_arms_keep_delimiters_hosts_and_endpoint_groups():
    source = "ac-C.[TBMB(4,4)[.!2(5,4)][.!3(6,4)]]-C.!2-C.!3-am"
    sequence = _sequence(source)
    peptide = Peptide.from_sequence(sequence)
    outer, first, second = peptide.layout.groups
    assert [(arm.host, arm.members, arm.parent) for arm in (first, second)] == [
        (outer.members[0], (), outer.id),
        (outer.members[0], (), outer.id),
    ]
    assert [(arm.opening, arm.closing) for arm in (first, second)] == [
        ("[", "]"),
        ("[", "]"),
    ]
    assert [(marker.label, marker.group) for marker in peptide.layout.markers] == [
        ("!2", first.id),
        ("!3", second.id),
        ("!2", None),
        ("!3", None),
    ]
    for notation in ("preserve", "percent", "bracket"):
        _assert_serialization(sequence, peptide, notation)


def test_new_layout_does_not_publish_stale_source_spans():
    peptide = Peptide.from_sequence(_sequence("ac-C.{!1(4,4)}-A-C.!1-am"))
    emitted = serialize(peptide, "bracket")
    assert emitted.layout.source is None
    assert all(group.span is None for group in emitted.layout.groups)
    assert all(marker.span is None for marker in emitted.layout.markers)


def test_provenance_cannot_claim_an_anchor_outside_its_occurrence():
    peptide = Peptide.from_sequence(_sequence("G"))
    wrong = replace(
        peptide.occurrences[0], provenance=AtomProvenance((0, 1), ((1, 2),), True)
    )
    with pytest.raises(ValueError, match="owned atom"):
        Peptide((wrong,), ())


def test_slot_identity_does_not_collapse_shared_anchor_atoms():
    peptide = Peptide.from_sequence(_sequence("G"))
    first, third = peptide.site(Endpoint(0, 1)), peptide.site(Endpoint(0, 3))
    assert first.anchor == third.anchor
    assert first.slot != third.slot
    assert peptide.connection_at(Endpoint(0, 1)) is None
    with pytest.raises(ValueError, match="more than once"):
        Peptide(
            peptide.occurrences,
            (
                Connection(Endpoint(0, 1), Endpoint(0, 2)),
                Connection(Endpoint(0, 3), Endpoint(0, 2)),
            ),
        )


def test_occupancy_reads_current_slots_and_connections_without_chemistry(monkeypatch):
    import pyPept.peptide as model

    sequence = Sequence("G-G")

    def unnecessary_enrichment(*args, **kwargs):
        raise AssertionError("Occupancy must not require chemical enrichment")

    monkeypatch.setattr(model, "attachment_sites", unnecessary_enrichment)
    assert Peptide.occupied_sites_from_sequence(sequence) == {
        Endpoint(0, 2),
        Endpoint(1, 1),
    }
    # R1 and R3 share an atom but have separate occupancy. Detached parser data
    # remains mutable, and inspection must validate its current contents.
    sequence.s_bonds[0][5] = 3
    assert Peptide.occupied_sites_from_sequence(sequence) == {
        Endpoint(0, 2),
        Endpoint(1, 3),
    }
    sequence.s_bonds.append(sequence.s_bonds[0][:])
    with pytest.raises(ValueError, match="more than once"):
        Peptide.occupied_sites_from_sequence(sequence)
    sequence.s_bonds.pop()
    sequence.s_bonds[0][5] = 99
    with pytest.raises(ValueError, match="missing attachment"):
        Peptide.occupied_sites_from_sequence(sequence)


def test_definition_reuse_is_local_to_each_projection(monkeypatch):
    import pyPept.peptide as model

    sequence = Sequence("G-G-G")
    discover = model.attachment_sites
    calls = []

    def counted(molecule, leaving_groups):
        calls.append((molecule, tuple(leaving_groups)))
        return discover(molecule, leaving_groups)

    monkeypatch.setattr(model, "attachment_sites", counted)
    first = Peptide.from_sequence(sequence)
    assert len(calls) == 1
    assert all(item.has_backbone for item in first.occurrences)
    assert all(
        {site.slot for site in item.sites} == {1, 2, 3} for item in first.occurrences
    )

    # The next projection must observe mutation of the same detached molecule.
    dummy = next(
        atom
        for atom in sequence.s_monomers[0]["m_romol"].GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 3
    )
    dummy.SetIsotope(7)
    second = Peptide.from_sequence(sequence)
    assert len(calls) == 2
    assert all(
        {site.slot for site in item.sites} == {1, 2, 7} for item in second.occurrences
    )
    assert all(
        {site.slot for site in item.sites} == {1, 2, 3} for item in first.occurrences
    )


def test_repeated_monomers_retain_occurrence_ids_and_input_atom_provenance():
    sequence = _sequence("A%G-G.!1(2,4)%K.!1(4,2)-A")
    peptide = Peptide.from_sequence(sequence)
    assembled = Molecule(sequence, depiction=None)
    atoms = assembled.get_residue_atom_map()
    occurrences = tuple(
        replace(item, provenance=AtomProvenance(tuple(atoms[item.id]), (), True))
        for item in peptide.occurrences
    )
    # This producer has no original source layout, as with a recognized partition.
    peptide = Peptide(occurrences, peptide.connections)
    percent = _assert_serialization(sequence, peptide, "percent")
    bracket = _assert_serialization(sequence, peptide, "bracket")
    assert percent.occurrence_order != bracket.occurrence_order
    for emission in (percent, bracket):
        assert set(emission.occurrence_order) == {item.id for item in occurrences}
        assigned = [
            peptide.occurrence(identity).provenance.source_atoms
            for identity in emission.occurrence_order
        ]
        assert sorted(atom for owned in assigned for atom in owned) == list(
            range(assembled.mol.GetNumAtoms())
        )


@pytest.mark.parametrize("notation", ["percent", "bracket"])
@pytest.mark.parametrize(
    "source",
    [
        "ac-Pra.!1(4,4)-A-A-AzK.!1(4,4)-am",
        "ac-S5.!1(4,4)-A-A-R8.!1(4,4)-am",
        "ac-K.[G(4,2).ac(1,2)]-C.!1(4,4)-am%ac-C.!1-G-am",
    ],
)
def test_reverse_import_uses_the_same_layout_with_exact_reaction_products(
    source, notation
):
    original = _molecule(_sequence(source))
    converted = convert_smiles(Chem.MolToSmiles(original), notation=notation)
    parsed = _sequence(converted.cabiln)
    assert converted.recognition_status == "complete"
    assert len(converted.assignments) == len(parsed.s_monomers)
    assert {item.residue_index: item.symbol for item in converted.assignments} == {
        index: item["m_abbr"] for index, item in enumerate(parsed.s_monomers)
    }
    assert sorted(
        atom for item in converted.assignments for atom in item.source_atoms
    ) == list(range(original.GetNumAtoms()))
    assert Chem.MolToSmiles(_molecule(parsed)) == Chem.MolToSmiles(original)
