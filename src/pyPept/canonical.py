"""Canonical encoding of one resolved monomer graph, without redecomposition.

The versioned convention uses RDKit to label a colored monomer/port graph. Its
colors are full definition identities and numbered slots, never source IDs or
source atom indices. The two writers share that order and preserve every
occurrence. A final parse supplies actual output indices and fresh layout spans.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

from rdkit import Chem, rdBase

from pyPept.peptide import (
    Connection, Endpoint, Peptide, Serialization,
)
from pyPept.source import Span


def canonical_convention():
    """Qualify exported text; RDKit upgrades can change canonical ordering."""
    return {
        "format": "cabiln-graph-v1",
        "labeling": "rdkit-colored-port-graph-v1",
        "rdkit": rdBase.rdkitVersion,
    }


@dataclass(frozen=True)
class _Component:
    order: tuple[int, ...]
    edges: tuple[Connection, ...]
    code: tuple


def _rank_component(identities, edges, nodes):
    graph = Chem.RWMol()
    symbols = []

    def vertex(color):
        atom = Chem.Atom(0)
        atom.SetNoImplicit(True)
        index = graph.AddAtom(atom)
        symbols.append(color)
        return index

    indices = {}
    for identity in identities:
        node = nodes[identity]
        color = json.dumps(
            (node.definition.key, node.has_backbone), separators=(",", ":")
        )
        indices[identity] = vertex("M:" + color)
    for edge in edges:
        a, b = edge.endpoints
        left = vertex(f"P:{a.slot}")
        right = vertex(f"P:{b.slot}")
        for start, end in (
            (indices[a.occurrence_id], left),
            (left, right),
            (right, indices[b.occurrence_id]),
        ):
            graph.AddBond(start, end, Chem.BondType.SINGLE)
    molecule = graph.GetMol()
    molecule.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(molecule)
    ranks = Chem.CanonicalRankAtomsInFragment(
        molecule,
        list(range(molecule.GetNumAtoms())),
        list(range(molecule.GetNumBonds())),
        symbols,
        breakTies=True,
        includeChirality=False,
        includeIsotopes=False,
        includeAtomMaps=False,
    )
    order = tuple(sorted(identities, key=lambda identity: ranks[indices[identity]]))
    positions = {identity: index for index, identity in enumerate(order)}
    code = (
        tuple((nodes[i].definition.key, nodes[i].has_backbone) for i in order),
        tuple(sorted(_edge_key(edge, positions) for edge in edges)),
    )
    return _Component(order, tuple(edges), code)


def _edge_key(edge, positions):
    return tuple(sorted(
        (positions[endpoint.occurrence_id], endpoint.slot)
        for endpoint in edge.endpoints
    ))


def _components(peptide, nodes):
    neighbors = {identity: [] for identity in nodes}
    for edge in peptide.connections:
        a, b = edge.endpoints
        neighbors[a.occurrence_id].append((b.occurrence_id, edge))
        neighbors[b.occurrence_id].append((a.occurrence_id, edge))
    remaining = set(nodes)
    result = []
    while remaining:
        stack = [next(iter(remaining))]
        members, edges = set(), set()
        while stack:
            identity = stack.pop()
            if identity in members:
                continue
            members.add(identity)
            for other, edge in neighbors[identity]:
                edges.add(edge)
                stack.append(other)
        remaining.difference_update(members)
        result.append(_rank_component(tuple(members), tuple(edges), nodes))
    return sorted(result, key=lambda component: component.code)


def _paths(component, nodes):
    positions = {identity: index for index, identity in enumerate(component.order)}
    forward, previous = {}, {}
    for edge in component.edges:
        a, b = edge.endpoints
        if (a.slot, b.slot) == (1, 2):
            a, b = b, a
        if (a.slot, b.slot) == (2, 1):
            forward[a.occurrence_id] = b.occurrence_id
            previous[b.occurrence_id] = a.occurrence_id
    remaining, paths = set(component.order), []
    starts = [identity for identity in component.order if identity not in previous]
    # Heads come first; remaining vertices belong to cycles. Their canonical
    # order fixes each cycle's cut without rescanning the remaining set.
    for current in starts + list(component.order):
        if current not in remaining:
            continue
        path = []
        while current in remaining:
            path.append(current)
            remaining.remove(current)
            current = forward.get(current)
        paths.append(tuple(path))
    return sorted(paths, key=lambda path: (
        -sum(nodes[i].has_backbone for i in path),
        -len(path), tuple(positions[i] for i in path),
    ))


def _tree(component, primary):
    positions = {identity: index for index, identity in enumerate(component.order)}
    neighbors = {identity: [] for identity in component.order}
    for edge in component.edges:
        a, b = edge.endpoints
        for own, other in ((a, b), (b, a)):
            neighbors[own.occurrence_id].append((own, other, edge))
    for entries in neighbors.values():
        entries.sort(key=lambda item: (
            item[0].slot, positions[item[1].occurrence_id], item[1].slot,
        ))
    discovered = set(primary)
    children, used = {}, set()
    for root in primary:
        stack = [(root, iter(neighbors[root]))]
        while stack:
            host, remaining = stack[-1]
            entry = next(remaining, None)
            if entry is None:
                stack.pop()
                continue
            own, other, edge = entry
            child = other.occurrence_id
            if child in discovered:
                continue
            discovered.add(child)
            used.add(edge)
            children.setdefault(host, []).append((child, own.slot, other.slot))
            stack.append((child, iter(neighbors[child])))
    if discovered != set(component.order):
        raise ValueError("Canonical tree did not visit every monomer occurrence")
    return children, used


def serialize_canonical(peptide, notation):
    """Emit the v1 convention and map actual parser indices back to occurrences."""
    if notation not in ("percent", "bracket"):
        raise ValueError("Canonical export requires 'percent' or 'bracket' notation")
    if not peptide.occurrences:
        raise ValueError("Cannot serialize an empty peptide")
    nodes = {item.id: item for item in peptide.occurrences}
    if any(node.definition is None for node in nodes.values()):
        raise ValueError(
            "Canonical export requires resolved monomer definitions; "
            "construct the peptide with Peptide.from_sequence"
        )
    plans, markers = [], {identity: [] for identity in nodes}
    next_tag = 1
    for component in _components(peptide, nodes):
        paths = _paths(component, nodes)
        children = {}
        if notation == "bracket":
            paths = paths[:1]
            children, consumed = _tree(component, paths[0])
        else:
            consumed = set()
        consumed.update(
            Connection(Endpoint(a, 2), Endpoint(b, 1))
            for path in paths for a, b in zip(path, path[1:])
        )
        positions = {
            identity: index for index, identity in enumerate(component.order)
        }
        for edge in sorted(
            (edge for edge in component.edges if edge not in consumed),
            key=lambda edge: _edge_key(edge, positions),
        ):
            a, b = edge.endpoints
            for own, other in ((a, b), (b, a)):
                markers[own.occurrence_id].append(
                    f".!{next_tag}({own.slot},{other.slot})"
                )
            next_tag += 1
        plans.extend((path, children) for path in paths)

    parts, spans = [], {}
    position = 0

    def append(text, identity=None):
        nonlocal position
        if identity is not None:
            spans[Span(position, position + len(text))] = identity
        parts.append(text)
        position += len(text)

    def write_node(root, children):
        stack = [("node", (root, None))]
        while stack:
            action, value = stack.pop()
            if action == "text":
                append(value)
                continue
            identity, pair = value
            append(nodes[identity].definition.token, identity)
            if pair is not None:
                append(f"({pair[0]},{pair[1]})")
            append("".join(markers[identity]))
            for child, own, other in reversed(children.get(identity, ())):
                stack.extend((
                    ("text", "]"), ("node", (child, (own, other))),
                    ("text", ".["),
                ))

    for index, (path, children) in enumerate(plans):
        if index:
            append("%")
        for offset, identity in enumerate(path):
            if offset:
                append("-")
            write_node(identity, children)
    text = "".join(parts)

    # Source spans identify the exact occurrence despite Sequence's separate
    # explicit/pendant collection order and arbitrary repeated symbols.
    from pyPept.sequence import Sequence

    parsed = Sequence(text, track_source=True, warning_sink=lambda _: None)
    order = tuple(spans[source.token] for source in parsed.s_sources)
    if len(order) != len(nodes) or set(order) != set(nodes):
        raise ValueError("Canonical serialization lost or repeated an occurrence")
    restored = Peptide.from_sequence(parsed, order)
    if set(restored.connections) != set(peptide.connections):
        raise ValueError("Canonical serialization changed monomer connections")
    for node in restored.occurrences:
        if node.definition.key != nodes[node.id].definition.key:
            raise ValueError("Canonical serialization changed a resolved definition")
    return Serialization(text, order, restored.layout)
