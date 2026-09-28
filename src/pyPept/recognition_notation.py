"""Turn a source-atom partition into CABILN without inferring a new partition.

Recognition decides atom ownership and attachment slots. This module preserves
that decision, fills explicitly unmatched regions, and chooses readable chains.
Every emitted interpretation still needs independent assembly verification.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from rdkit import Chem

from pyPept.peptide import (
    AtomProvenance,
    AttachmentSite,
    Connection,
    Endpoint,
    MonomerOccurrence,
    Peptide,
    serialize,
)

if TYPE_CHECKING:
    from pyPept.recognition import Candidate


@dataclass(frozen=True)
class MonomerAssignment:
    """An occurrence in the returned CABILN and its original source atoms.

    ``residue_index`` is Sequence's zero-based occurrence index. It is explicit
    because unresolved components can have occurrences without assignments.
    """

    symbol: str
    source_atoms: tuple[int, ...]
    attachments: tuple[tuple[int, int], ...]
    recognized: bool
    residue_index: int


@dataclass(frozen=True)
class Interpretation:
    cabiln: str
    assignments: tuple[MonomerAssignment, ...]
    details: tuple[tuple[str, float, int], ...]
    unknown_pieces: int
    recognized_atoms: int
    backbone_connections: int


class UnknownRegionError(ValueError):
    """An unknown region is impossible regardless of boundary slot choices."""

    def __init__(self, atoms, message):
        super().__init__(message)
        self.atoms = frozenset(atoms)


def _carbonyl(molecule, index):
    atom = molecule.GetAtomWithIdx(index)
    return atom.GetAtomicNum() == 6 and any(
        neighbor.GetAtomicNum() == 8
        and molecule.GetBondBetweenAtoms(index, neighbor.GetIdx()).GetBondType()
        == Chem.BondType.DOUBLE
        for neighbor in atom.GetNeighbors()
    )


def _regions(molecule, atoms):
    remaining = set(atoms)
    while remaining:
        pending = [min(remaining)]
        region = set()
        while pending:
            index = pending.pop()
            if index not in remaining:
                continue
            remaining.remove(index)
            region.add(index)
            pending.extend(
                atom.GetIdx()
                for atom in molecule.GetAtomWithIdx(index).GetNeighbors()
                if atom.GetIdx() in remaining
            )
        yield frozenset(region)


def _unknown_candidate(molecule, atoms, known_ports):
    """Keep one unmatched region and the source bonds at its known boundaries.

    R1/R2 at a boundary follow the neighbouring monomer's actual R2/R1. Free
    termini are taken from one unambiguous amine/carboxyl pair. Other boundaries
    are explicit sidechain ports. No missing sidechain is replaced with a name.
    """
    from pyPept.recognition import Candidate

    # Every supported local synthetic template needs a carbonyl R2 anchor.
    # This condition is independent of neighboring aliases and slot choices.
    # The search must separately prove that no known template can cover this
    # region before using the failure to reject other ownership partitions.
    if not any(_carbonyl(molecule, index) for index in atoms):
        raise UnknownRegionError(
            atoms, "An unmatched region has no carbonyl for an editable R2 site"
        )

    boundaries = []
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a in atoms) == (b in atoms):
            continue
        inside, outside = (a, b) if a in atoms else (b, a)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("An unmatched region has an unsupported boundary bond")
        partner_slot = known_ports.get((outside, inside))
        if partner_slot is None:
            raise ValueError("An unmatched region lacks a known attachment boundary")
        boundaries.append((inside, outside, partner_slot, bond.GetIdx()))
    boundaries.sort()

    # Backbone slots have chemical as well as positional meaning. A known
    # sidechain attachment alone does not establish an unknown residue boundary.
    incoming = [
        edge
        for edge in boundaries
        if edge[2] == 2 and molecule.GetAtomWithIdx(edge[0]).GetAtomicNum() in (7, 8)
    ]
    outgoing = [
        edge for edge in boundaries if edge[2] == 1 and _carbonyl(molecule, edge[0])
    ]
    if len(incoming) > 1 or len(outgoing) > 1:
        raise ValueError("Unmatched region has ambiguous backbone attachment sites")
    n_atom = incoming[0][0] if incoming else None
    c_atom = outgoing[0][0] if outgoing else None

    free_n = [
        i
        for i in atoms
        if molecule.GetAtomWithIdx(i).GetAtomicNum() == 7
        and molecule.GetAtomWithIdx(i).GetTotalNumHs() > 0
        and not molecule.GetAtomWithIdx(i).GetIsAromatic()
        and all(
            bond.GetBondType() == Chem.BondType.SINGLE
            for bond in molecule.GetAtomWithIdx(i).GetBonds()
        )
        and not any(
            _carbonyl(molecule, a.GetIdx())
            for a in molecule.GetAtomWithIdx(i).GetNeighbors()
        )
    ]
    acid_o = [
        i
        for i in atoms
        if molecule.GetAtomWithIdx(i).GetAtomicNum() == 8
        and molecule.GetAtomWithIdx(i).GetDegree() == 1
        and molecule.GetAtomWithIdx(i).GetTotalNumHs() == 1
        and _carbonyl(molecule, molecule.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx())
    ]
    acid_by_c = {
        molecule.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx(): i for i in acid_o
    }
    # The nearest free termini are the same topology convention used by monomer
    # ingestion. Require a unique minimum; a tie is not a discovered backbone.
    if n_atom is None or c_atom is None:
        pairs = []
        for n in [n_atom] if n_atom is not None else free_n:
            for c in [c_atom] if c_atom is not None else acid_by_c:
                path = Chem.GetShortestPath(molecule, n, c)
                if path and set(path) <= atoms:
                    pairs.append((len(path), n, c))
        pairs.sort()
        if not pairs or (len(pairs) > 1 and pairs[0][0] == pairs[1][0]):
            raise ValueError(
                "Cannot establish editable termini for an unmatched region"
            )
        _, n_atom, c_atom = pairs[0]

    edge_slots = {}
    next_side = 4
    for inside, outside, partner_slot, bond_index in boundaries:
        if incoming and (inside, outside) == incoming[0][:2]:
            slot = 1
        elif outgoing and (inside, outside) == outgoing[0][:2]:
            slot = 2
        else:
            slot, next_side = next_side, next_side + 1
        edge_slots[bond_index] = slot

    # FragmentOnBonds retains neighbour order and stereo references. Only the
    # original inter-region bonds are cut; interior source chemistry is kept.
    cuts = list(edge_slots)
    opened = (
        Chem.FragmentOnBonds(
            molecule, cuts, dummyLabels=[(edge_slots[i], edge_slots[i]) for i in cuts]
        )
        if cuts
        else Chem.Mol(molecule)
    )
    mappings = []
    fragments = Chem.GetMolFrags(opened, asMols=True, fragsMolAtomMapping=mappings)
    original = set(range(molecule.GetNumAtoms()))
    fragment, mapping = next(
        (fragment, mapping)
        for fragment, mapping in zip(fragments, mappings)
        if set(mapping) & original == set(atoms)
    )
    rw = Chem.RWMol(fragment)
    if not outgoing:
        oxygen = acid_by_c[c_atom]
        dummy = Chem.Atom(0)
        dummy.SetIsotope(2)
        rw.ReplaceAtom(mapping.index(oxygen), dummy)
    if not incoming:
        nitrogen = rw.GetAtomWithIdx(mapping.index(n_atom))
        if nitrogen.GetNumExplicitHs():
            nitrogen.SetNumExplicitHs(nitrogen.GetNumExplicitHs() - 1)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(1)
        index = rw.AddAtom(dummy)
        rw.AddBond(mapping.index(n_atom), index, Chem.BondType.SINGLE)
    template = rw.GetMol()
    Chem.SanitizeMol(template)
    token = f"<{Chem.MolToSmiles(template)}>"
    anchors = [(n_atom, 1), (c_atom, 2)]
    ports = []
    for inside, outside, _, bond_index in boundaries:
        slot = edge_slots[bond_index]
        ports.append((inside, outside, slot))
        if (inside, slot) not in anchors:
            anchors.append((inside, slot))
    return Candidate(
        symbol=token,
        atoms=atoms,
        ports=tuple(ports),
        anchors=tuple(anchors),
        has_backbone=True,
        kind="unknown",
        template_size=len(atoms),
        attachment_types=(),
    )


def _connect(nodes, extra_edges):
    owner = {a: i for i, node in enumerate(nodes) for a in node.atoms}
    if len(owner) != sum(len(node.atoms) for node in nodes):
        raise ValueError("A decomposition assigns a source atom more than once")
    ports = {
        (i, a, b): slot for i, node in enumerate(nodes) for a, b, slot in node.ports
    }
    edges = set()
    for (i, a, b), slot in ports.items():
        j = owner[b]
        partner = ports.get((j, b, a))
        if partner is None:
            raise ValueError("A decomposition connection lacks a reciprocal slot")
        edge = tuple(sorted(((i, slot), (j, partner))))
        edges.add(edge)
    for a, slot_a, b, slot_b in extra_edges:
        edges.add(tuple(sorted(((owner[a], slot_a), (owner[b], slot_b)))))
    return tuple(
        Connection(Endpoint(*left), Endpoint(*right)) for left, right in sorted(edges)
    )


def emit_interpretation(
    molecule,
    candidates: tuple[Candidate, ...],
    *,
    extra_edges=(),
    atom_origins=None,
    notation="percent",
):
    """Preserve a candidate cover, including unknown regions, in readable notation."""
    nodes = list(candidates)
    owned = set().union(*(node.atoms for node in nodes)) if nodes else set()
    all_atoms = set(range(molecule.GetNumAtoms()))
    if not owned <= all_atoms:
        raise ValueError("A decomposition references an atom outside its source")
    known_ports = {(a, b): slot for node in nodes for a, b, slot in node.ports}
    for atoms in _regions(molecule, all_atoms - owned):
        nodes.append(_unknown_candidate(molecule, atoms, known_ports))
    edges = _connect(nodes, extra_edges)
    origins = (
        tuple(range(molecule.GetNumAtoms())) if atom_origins is None else atom_origins
    )
    occurrences = tuple(
        MonomerOccurrence(
            id=i,
            symbol=node.symbol,
            sites=tuple(
                AttachmentSite(slot, dict(node.attachment_types).get(slot, ""))
                for slot in sorted({slot for _, slot in node.anchors})
            ),
            has_backbone=node.has_backbone,
            kind=node.kind,
            size=len(node.atoms),
            order_key=tuple(sorted(node.atoms)),
            provenance=AtomProvenance(
                source_atoms=tuple(
                    sorted(origins[a] for a in node.atoms if origins[a] is not None)
                ),
                attachments=tuple(
                    sorted(
                        (slot, origins[a])
                        for a, slot in node.anchors
                        if origins[a] is not None
                    )
                ),
                recognized=node.kind != "unknown",
            ),
        )
        for i, node in enumerate(nodes)
    )
    peptide = Peptide(occurrences, edges)
    emitted = serialize(peptide, notation=notation)
    assignments = tuple(
        MonomerAssignment(
            symbol=nodes[i].symbol,
            source_atoms=occurrences[i].provenance.source_atoms,
            attachments=occurrences[i].provenance.attachments,
            recognized=occurrences[i].provenance.recognized,
            residue_index=residue_index,
        )
        for residue_index, i in enumerate(emitted.occurrence_order)
    )
    # Keep the legacy residue-details view on the primary backbone. The complete
    # partition, including branches and caps, is recorded in assignments.
    primary = emitted.layout.segments[0].roots
    residue_indices = [i for i in primary if nodes[i].has_backbone]
    if not residue_indices:
        residue_indices = list(primary)
    details = tuple(
        (
            nodes[i].symbol,
            len(nodes[i].atoms) if nodes[i].kind != "unknown" else 0,
            len(nodes[i].atoms),
        )
        for i in residue_indices
    )
    return Interpretation(
        cabiln=emitted.text,
        assignments=assignments,
        details=details,
        unknown_pieces=sum(node.kind == "unknown" for node in nodes),
        recognized_atoms=sum(len(a.source_atoms) for a in assignments if a.recognized),
        backbone_connections=sum(
            {edge.left.slot, edge.right.slot} == {1, 2}
            and nodes[edge.left.occurrence_id].has_backbone
            and nodes[edge.right.occurrence_id].has_backbone
            for edge in edges
        ),
    )


def opaque_encodings(fragment, ring_tag=1, max_breaks=16):
    """Yield full-structure synthetic candidates without guessing residues.

    Try replacing the termini of an amino acid first. Otherwise open an amide:
    a ring becomes one synthetic residue with a terminal crosslink, while an
    acyclic molecule becomes two synthetic caps joined by the original amide.
    Each candidate must still be checked by assembly against the source.
    """
    from rdkit import Chem

    def carbonyl(atom):
        return atom.GetAtomicNum() == 6 and any(
            nb.GetAtomicNum() == 8
            and fragment.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondType()
            == Chem.BondType.DOUBLE
            for nb in atom.GetNeighbors()
        )

    free_nitrogens = [
        atom.GetIdx()
        for atom in fragment.GetAtoms()
        if atom.GetAtomicNum() == 7
        and atom.GetTotalNumHs() >= 1
        and not any(carbonyl(nb) for nb in atom.GetNeighbors())
    ]
    acid_oxygens = [
        atom.GetIdx()
        for atom in fragment.GetAtoms()
        if atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and carbonyl(atom.GetNeighbors()[0])
    ]
    if free_nitrogens and acid_oxygens:
        # Replace the acid OH; adding another substituent to its carbonyl C
        # would produce pentavalent carbon (the former fallback did exactly that).
        rw = Chem.RWMol(fragment)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(2)
        rw.ReplaceAtom(acid_oxygens[0], dummy)
        n_idx = free_nitrogens[0]
        nitrogen = rw.GetAtomWithIdx(n_idx)
        if nitrogen.GetNumExplicitHs():
            nitrogen.SetNumExplicitHs(nitrogen.GetNumExplicitHs() - 1)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(1)
        d_idx = rw.AddAtom(dummy)
        rw.AddBond(n_idx, d_idx, Chem.BondType.SINGLE)
        try:
            molecule = rw.GetMol()
            Chem.SanitizeMol(molecule)
            yield f"<{Chem.MolToSmiles(molecule)}>", False
        except (ValueError, RuntimeError):
            pass

    n_tried = 0
    for bond in fragment.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 7 and carbonyl(b):
            nitrogen, carbon = a, b
        elif b.GetAtomicNum() == 7 and carbonyl(a):
            nitrogen, carbon = b, a
        else:
            continue
        if n_tried >= max_breaks:
            break
        n_tried += 1
        rw = Chem.RWMol(fragment)
        rw.RemoveBond(nitrogen.GetIdx(), carbon.GetIdx())
        for slot, atom in ((1, nitrogen), (2, carbon)):
            dummy = Chem.Atom(0)
            dummy.SetIsotope(slot)
            d_idx = rw.AddAtom(dummy)
            rw.AddBond(atom.GetIdx(), d_idx, Chem.BondType.SINGLE)
        try:
            opened = rw.GetMol()
            Chem.SanitizeMol(opened)
            parts = Chem.GetMolFrags(opened, asMols=True)
            if len(parts) == 1:
                token = Chem.MolToSmiles(parts[0])
                yield f"!{ring_tag}-<{token}>-!{ring_tag}", True
            elif len(parts) == 2:
                caps = {}
                for part in parts:
                    slot = next(
                        atom.GetIsotope()
                        for atom in part.GetAtoms()
                        if atom.GetAtomicNum() == 0
                    )
                    caps[slot] = Chem.MolToSmiles(part)
                yield f"<{caps[2]}>_-_<{caps[1]}>", False
        except (ValueError, RuntimeError):
            continue
