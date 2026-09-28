"""Resolved monomer occurrences, numbered connections, and notation layout.

Occurrence IDs survive formatting; output residue indices do not. Sequence and
recognition are adapters into this model. Chemistry remains in Molecule and its
reaction definitions. Source grouping is recorded once, never inferred from
which chemical component happens to be the first chain.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace

from pyPept.attachments import attachment_sites
from pyPept.source import Occurrence as SourceOccurrence
from pyPept.source import Span


@dataclass(frozen=True, order=True)
class Endpoint:
    occurrence_id: int
    slot: int


@dataclass(frozen=True)
class AttachmentSite:
    slot: int
    chem_type: str = ""
    anchor: int | None = None  # Template-local atom index, never a source atom ID.


@dataclass(frozen=True)
class AtomProvenance:
    source_atoms: tuple[int, ...]
    attachments: tuple[tuple[int, int], ...]  # (slot, original input atom)
    recognized: bool


@dataclass(frozen=True)
class MonomerOccurrence:
    id: int
    symbol: str
    sites: tuple[AttachmentSite, ...]
    has_backbone: bool = False
    kind: str = "aa"
    size: int = 0
    token: str | None = None
    source: SourceOccurrence | None = None
    provenance: AtomProvenance | None = None
    order_key: tuple[int, ...] = ()


@dataclass(frozen=True, order=True)
class Connection:
    left: Endpoint
    right: Endpoint
    label: str | None = field(default=None, compare=False)

    def __post_init__(self):
        if self.right < self.left:
            old_left = self.left
            object.__setattr__(self, "left", self.right)
            object.__setattr__(self, "right", old_left)

    @property
    def endpoints(self):
        return self.left, self.right


@dataclass(frozen=True)
class LayoutSegment:
    id: int
    roots: tuple[int, ...]
    members: tuple[int, ...]


@dataclass(frozen=True)
class LayoutGroup:
    id: int
    host: int
    members: tuple[int, ...]
    kind: str  # inline, bracket, or arm
    opening: str
    closing: str
    protected: bool
    span: Span | None = None
    parent: int | None = None


@dataclass(frozen=True)
class LayoutMarker:
    endpoint: Endpoint
    label: str
    kind: str  # inline, bracket, or terminal
    span: Span | None = None
    group: int | None = None


@dataclass(frozen=True)
class NotationLayout:
    """Presentation facts; ``source`` exists only with valid original locations.

    A serializer-generated layout has no original source spans. Its text and
    output occurrence order are provided by Serialization instead.
    """

    source: str | None
    segments: tuple[LayoutSegment, ...]
    groups: tuple[LayoutGroup, ...] = ()
    markers: tuple[LayoutMarker, ...] = ()


@dataclass(frozen=True)
class Peptide:
    occurrences: tuple[MonomerOccurrence, ...]
    connections: tuple[Connection, ...]
    layout: NotationLayout | None = None

    def __post_init__(self):
        identities = {item.id for item in self.occurrences}
        if len(identities) != len(self.occurrences):
            raise ValueError("Monomer occurrence IDs must be unique")
        sites = set()
        atoms = set()
        for item in self.occurrences:
            for site in item.sites:
                endpoint = Endpoint(item.id, site.slot)
                if site.slot < 1 or endpoint in sites:
                    raise ValueError(
                        "Attachment sites need unique positive slot numbers"
                    )
                sites.add(endpoint)
            if item.provenance is not None:
                owned = set(item.provenance.source_atoms)
                if len(owned) != len(item.provenance.source_atoms) or atoms & owned:
                    raise ValueError(
                        "A source atom belongs to more than one occurrence"
                    )
                atoms.update(owned)
                if any(
                    atom not in owned or Endpoint(item.id, slot) not in sites
                    for slot, atom in item.provenance.attachments
                ):
                    raise ValueError(
                        "Source attachment must name an owned atom and an existing slot"
                    )
        used = set()
        labels = set()
        for connection in self.connections:
            if connection.label is not None:
                if (
                    not re.fullmatch(r"!\w+", connection.label)
                    or connection.label in labels
                ):
                    raise ValueError(
                        "Crosslink labels must be valid and unique per connection"
                    )
                labels.add(connection.label)
            for endpoint in connection.endpoints:
                if endpoint not in sites:
                    raise ValueError(
                        f"Connection refers to a missing attachment: {endpoint}"
                    )
                if endpoint in used:
                    raise ValueError(
                        f"Attachment site is used more than once: {endpoint}"
                    )
                used.add(endpoint)
        if self.layout is not None:
            group_ids = {group.id for group in self.layout.groups}
            if len(group_ids) != len(self.layout.groups):
                raise ValueError("Layout group IDs must be unique")
            for segment in self.layout.segments:
                if (
                    not set(segment.roots) <= identities
                    or not set(segment.members) <= identities
                ):
                    raise ValueError("Layout segment refers to a missing occurrence")
            for group in self.layout.groups:
                if group.host not in identities or not set(group.members) <= identities:
                    raise ValueError("Layout group refers to a missing occurrence")
                if group.parent is not None and (
                    group.parent not in group_ids or group.parent == group.id
                ):
                    raise ValueError("Layout group refers to an invalid parent")

    def occurrence(self, identity):
        for item in self.occurrences:
            if item.id == identity:
                return item
        raise ValueError(f"Monomer occurrence {identity} does not exist")

    def site(self, endpoint):
        for site in self.occurrence(endpoint.occurrence_id).sites:
            if site.slot == endpoint.slot:
                return site
        raise ValueError(f"Attachment site {endpoint} does not exist")

    def connection_at(self, endpoint):
        self.site(endpoint)
        return next(
            (edge for edge in self.connections if endpoint in edge.endpoints), None
        )

    @classmethod
    def from_sequence(cls, sequence, occurrence_ids=None):
        """Adapt a parsed Sequence; source layout requires ``track_source=True``.

        ``occurrence_ids`` optionally maps Sequence indices to retained IDs after
        an edit. No symbol matching or molecule isomorphism establishes identity.
        """
        ids = (
            tuple(range(len(sequence.s_monomers)))
            if occurrence_ids is None
            else tuple(occurrence_ids)
        )
        if len(ids) != len(sequence.s_monomers):
            raise ValueError("Expected one occurrence ID per Sequence monomer")
        tracked = len(sequence.s_sources) == len(ids)
        source = str(sequence.s_inputbiln) if tracked else None
        occurrences = []
        for index, monomer in enumerate(sequence.s_monomers):
            mol = monomer["m_romol"]
            anchors = {
                atom.GetIsotope(): atom.GetNeighbors()[0].GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 0
            }
            sites = tuple(
                AttachmentSite(site["slot"], site["chem_type"], anchors[site["slot"]])
                for site in attachment_sites(mol, monomer["m_Rgroups"])
            )
            types = {site.slot: site.chem_type for site in sites}
            location = sequence.s_sources[index] if tracked else None
            occurrences.append(
                MonomerOccurrence(
                    ids[index],
                    monomer["m_abbr"],
                    sites,
                    types.get(1) in ("backbone_n", "backbone_o")
                    and types.get(2)
                    in ("backbone_c", "backbone_c_red", "sp3_c_anchor"),
                    monomer["m_type"],
                    mol.GetNumHeavyAtoms(),
                    (
                        source[location.token.start : location.token.end]
                        if tracked
                        else None
                    ),
                    location,
                )
            )
        connections = tuple(
            Connection(Endpoint(ids[b[0]], b[4]), Endpoint(ids[b[2]], b[5]))
            for b in sequence.s_bonds
        )
        if not tracked:
            return cls(tuple(occurrences), connections)
        layout = _source_layout(
            source, tuple(occurrences), sequence.s_biln.tracker.markers
        )
        labels = {marker.endpoint: marker.label for marker in layout.markers}
        connections = tuple(
            Connection(
                edge.left, edge.right, labels.get(edge.left, labels.get(edge.right))
            )
            for edge in connections
        )
        return cls(tuple(occurrences), connections, layout)


def _source_layout(source, occurrences, markers):
    """Project parser-recorded groups and owners, including all source segments."""
    source_groups = {}
    for item in occurrences:
        location = item.source
        if location.bracket is not None:
            source_groups[location.bracket] = "bracket"
        if location.arm is not None:
            source_groups[location.arm] = "arm"
        if location.kind == "inline":
            source_groups[location.entry] = "inline"
    for marker in markers:
        if marker.bracket is not None:
            source_groups[marker.bracket] = "bracket"
        if marker.arm is not None:
            source_groups[marker.arm] = "arm"
    spans = sorted(source_groups, key=lambda span: (span.start, -span.end))
    members = {
        span: tuple(
            item.id
            for item in occurrences
            if span.start <= item.source.token.start
            and item.source.token.end <= span.end
        )
        for span in spans
    }
    parents = {}
    for span in spans:
        enclosing = [
            other
            for other in spans
            if other != span and other.start <= span.start and span.end <= other.end
        ]
        parents[span] = (
            min(enclosing, key=lambda other: other.end - other.start)
            if enclosing
            else None
        )
    groups = []
    for identity, span in enumerate(spans):
        parent = parents[span]
        if parent is None:
            candidates = [
                item
                for item in occurrences
                if item.id not in members[span]
                and item.source.entry.start <= span.start < item.source.entry.end
            ]
            candidates.sort(
                key=lambda item: item.source.entry.end - item.source.entry.start
            )
        else:
            sibling_members = {
                member
                for other in spans
                if parents[other] == parent
                for member in members[other]
            }
            candidates = [
                item
                for item in occurrences
                if item.id in members[parent]
                and item.id not in sibling_members
                and item.source.token.end <= span.start
            ]
            candidates.sort(key=lambda item: -item.source.token.end)
        if not candidates:
            raise ValueError("Source group has no recorded host occurrence")
        raw = source[span.start : span.end].lstrip(".")
        kind = source_groups[span]
        opening = "" if kind == "inline" else raw[0]
        if kind == "arm" and opening not in ("[", "{"):
            # Legacy [[anchor.arm].arm] lowers to sibling arms whose original
            # spans cover entries, not bracket characters. The semantic layout
            # uses square arms while preserve-mode retains the literal source.
            opening = "["
        groups.append(
            LayoutGroup(
                identity,
                candidates[0].id,
                members[span],
                kind,
                opening,
                {"[": "]", "{": "}"}.get(opening, ""),
                opening == "{",
                span,
                spans.index(parent) if parent is not None else None,
            )
        )
    owner_ids = {
        span: item.id
        for item in occurrences
        for span in (item.source.token, item.source.entry)
    }
    resolved_markers = []
    for marker in markers:
        match = re.search(r"!\w+", source[marker.span.start : marker.span.end])
        if match is None or marker.owner not in owner_ids:
            raise ValueError("Crosslink marker lost its source label or owner")
        containing = [
            span
            for span in spans
            if span.start <= marker.span.start and marker.span.end <= span.end
        ]
        group = (
            spans.index(min(containing, key=lambda span: span.end - span.start))
            if containing
            else None
        )
        resolved_markers.append(
            LayoutMarker(
                Endpoint(owner_ids[marker.owner], marker.slot),
                match.group(),
                marker.kind,
                marker.span,
                group,
            )
        )
    segment_ids = sorted({item.source.segment for item in occurrences})
    segments = tuple(
        LayoutSegment(
            segment,
            tuple(
                item.id
                for item in occurrences
                if item.source.segment == segment and item.source.kind == "explicit"
            ),
            tuple(item.id for item in occurrences if item.source.segment == segment),
        )
        for segment in segment_ids
    )
    return NotationLayout(source, segments, tuple(groups), tuple(resolved_markers))


@dataclass(frozen=True)
class Serialization:
    text: str
    occurrence_order: tuple[int, ...]
    layout: NotationLayout


def _chains(peptide, identities, ignored):
    """Select R2-to-R1 paths only after occurrence connections are fixed."""
    forward, previous = {}, {}
    for edge in peptide.connections:
        if edge in ignored:
            continue
        a, b = edge.endpoints
        if a.occurrence_id not in identities or b.occurrence_id not in identities:
            continue
        if (a.slot, b.slot) == (1, 2):
            a, b = b, a
        if (a.slot, b.slot) == (2, 1):
            forward[a.occurrence_id] = b.occurrence_id
            previous[b.occurrence_id] = a.occurrence_id

    def key(identity):
        node = peptide.occurrence(identity)
        return (
            -node.size,
            len(node.symbol),
            node.symbol,
            node.order_key or (identity,),
        )

    remaining, chains = set(identities), []
    while remaining:
        starts = [identity for identity in remaining if identity not in previous]
        start = min(starts or remaining, key=key)
        chain, at = [], start
        while at in remaining:
            chain.append(at)
            remaining.remove(at)
            at = forward.get(at)
        chains.append((tuple(chain), at == start))
    chains.sort(
        key=lambda item: (
            -sum(peptide.occurrence(i).has_backbone for i in item[0]),
            -len(item[0]),
            tuple(key(i) for i in item[0]),
        )
    )
    return chains


def _fixed_source_groups(peptide, notation):
    layout = peptide.layout
    if layout is None or layout.source is None:
        return ()
    result = []
    for group in layout.groups:
        if group.parent is not None:
            continue
        nested = any(child.parent == group.id for child in layout.groups)
        members = set(group.members)
        backbone = sum(
            {edge.left.slot, edge.right.slot} == {1, 2}
            and all(endpoint.occurrence_id in members for endpoint in edge.endpoints)
            for edge in peptide.connections
        )
        if (
            notation == "bracket"
            or group.protected
            or nested
            or group.kind == "inline"
            or len(members) < 2
            or backbone != len(members) - 1
        ):
            result.append(group)
    return tuple(result)


def serialize(peptide: Peptide, notation="percent") -> Serialization:
    """Format a resolved peptide and return exact output-index to occurrence IDs.

    ``preserve`` returns tracked source verbatim. ``percent`` and ``bracket``
    choose readable layouts without changing connections. Protected groups,
    nested source arms and unsupported bracket tokens retain valid source forms.
    Bracket formatting may therefore retain explicit segments. It is a layout
    preference, not permission to alter chemistry or protected syntax.
    """
    if notation not in ("preserve", "percent", "bracket"):
        raise ValueError(f"Unknown notation layout: {notation!r}")
    if not peptide.occurrences:
        raise ValueError("Cannot serialize an empty peptide")
    if notation == "preserve":
        if peptide.layout is None or peptide.layout.source is None:
            raise ValueError("This peptide has no original source to preserve")
        return Serialization(
            peptide.layout.source,
            tuple(item.id for item in peptide.occurrences),
            peptide.layout,
        )

    from pyPept.sequence import _BRACKET_ENTRY_RE

    nodes = {item.id: item for item in peptide.occurrences}
    fixed = _fixed_source_groups(peptide, notation)
    fixed_members = {member for group in fixed for member in group.members}
    fixed_at = {}
    represented = set()
    for group in fixed:
        fixed_at.setdefault(group.host, []).append(group)
        members = set(group.members)
        marked = {
            marker.endpoint
            for marker in peptide.layout.markers
            if marker.span is not None
            and group.span.start <= marker.span.start
            and marker.span.end <= group.span.end
        }
        represented.update(marked)
        for edge in peptide.connections:
            inside = [
                endpoint
                for endpoint in edge.endpoints
                if endpoint.occurrence_id in members
            ]
            if len(inside) == 2:
                represented.update(edge.endpoints)
            elif len(inside) == 1:
                outside = next(
                    endpoint for endpoint in edge.endpoints if endpoint != inside[0]
                )
                if outside.occurrence_id == group.host and inside[0] not in marked:
                    represented.update(edge.endpoints)
                elif inside[0] not in marked:
                    raise ValueError(
                        "Preserved source group lacks its external connection marker"
                    )
    consumed = {
        edge for edge in peptide.connections if set(edge.endpoints) <= represented
    }
    chains = _chains(peptide, set(nodes) - fixed_members, consumed)
    chain_of = {
        identity: ci for ci, (chain, _) in enumerate(chains) for identity in chain
    }
    attachments = {}
    if notation == "bracket":
        for ci, (chain, cyclic) in enumerate(chains):
            if cyclic or any(identity in fixed_at for identity in chain):
                continue
            if not all(
                _BRACKET_ENTRY_RE.fullmatch(
                    f".{nodes[i].token or nodes[i].symbol}(1,1)"
                )
                for i in chain
            ):
                continue
            options = []
            for edge in peptide.connections:
                for own, host in (edge.endpoints, tuple(reversed(edge.endpoints))):
                    parent = chain_of.get(host.occurrence_id)
                    if (
                        own.occurrence_id in chain
                        and parent is not None
                        and parent < ci
                        and parent not in attachments
                        and edge not in consumed
                    ):
                        options.append(
                            (
                                parent,
                                chains[parent][0].index(host.occurrence_id),
                                chain.index(own.occurrence_id),
                                edge,
                                host,
                                own,
                            )
                        )
            if options:
                _, _, _, edge, host, own = min(options)
                attachments[ci] = host, own, edge

    linear_edges = {
        Connection(Endpoint(a, 2), Endpoint(b, 1))
        for chain, _ in chains
        for a, b in zip(chain, chain[1:])
    }
    consumed.update(linear_edges)
    consumed.update(edge for _, _, edge in attachments.values())
    reserved = {edge.label for edge in peptide.connections if edge.label}
    labels = {}
    next_label = 1

    def label(edge):
        nonlocal next_label
        if edge in labels:
            return labels[edge]
        if edge.label:
            value = edge.label
        else:
            while f"!{next_label}" in reserved:
                next_label += 1
            value = f"!{next_label}"
            reserved.add(value)
            next_label += 1
        labels[edge] = value
        return value

    cyclic_tags = {}
    for ci, (chain, cyclic) in enumerate(chains):
        if cyclic:
            closing = Connection(Endpoint(chain[-1], 2), Endpoint(chain[0], 1))
            edge = next(item for item in peptide.connections if item == closing)
            consumed.add(edge)
            cyclic_tags[ci] = label(edge)
    suffixes = {identity: [] for identity in nodes}
    generated_markers = []
    for edge in sorted(peptide.connections):
        if edge in consumed:
            continue
        tag = label(edge)
        for own, partner in (edge.endpoints, tuple(reversed(edge.endpoints))):
            if own in represented:
                continue
            if own.occurrence_id in fixed_members:
                raise ValueError(
                    "A fixed source group cannot receive a generated endpoint"
                )
            suffixes[own.occurrence_id].append(f".{tag}({own.slot},{partner.slot})")
            generated_markers.append(LayoutMarker(own, tag, "inline"))

    groups = []
    if fixed:
        for group in peptide.layout.groups:
            if any(
                root.span.start <= group.span.start and group.span.end <= root.span.end
                for root in fixed
            ):
                groups.append(
                    LayoutGroup(
                        group.id,
                        group.host,
                        group.members,
                        group.kind,
                        group.opening,
                        group.closing,
                        group.protected,
                        None,
                        group.parent,
                    )
                )
    next_group = max((group.id for group in groups), default=-1) + 1
    branches_at = {}
    for ci, (host, own, edge) in attachments.items():
        branches_at.setdefault(host.occurrence_id, []).append((ci, host, own))
    explicit_order, pendant_order = [], []

    def entry(identity, previous_slot, own_slot):
        return (
            f"{nodes[identity].token or nodes[identity].symbol}({previous_slot},{own_slot})"
            + "".join(suffixes[identity])
        )

    def branch_text(ci, host, own):
        nonlocal next_group
        chain = chains[ci][0]
        at = chain.index(own.occurrence_id)
        after, before = chain[at + 1 :], tuple(reversed(chain[:at]))
        order = (own.occurrence_id,) + after + before
        pendant_order.extend(order)
        identity = next_group
        next_group += 1
        groups.append(
            LayoutGroup(identity, host.occurrence_id, order, "bracket", "[", "]", False)
        )
        text = entry(own.occurrence_id, host.slot, own.slot)
        if after and before:
            text += "[." + ".".join(entry(i, 2, 1) for i in after) + "]"
            groups.append(
                LayoutGroup(
                    next_group,
                    own.occurrence_id,
                    after,
                    "arm",
                    "[",
                    "]",
                    False,
                    parent=identity,
                )
            )
            next_group += 1
        elif after:
            text += "." + ".".join(entry(i, 2, 1) for i in after)
        if before:
            text += "." + ".".join(entry(i, 1, 2) for i in before)
        return ".[" + text + "]"

    def token(identity):
        explicit_order.append(identity)
        text = nodes[identity].token or nodes[identity].symbol
        for group in sorted(
            fixed_at.get(identity, ()), key=lambda item: item.span.start
        ):
            text += peptide.layout.source[group.span.start : group.span.end]
            pendant_order.extend(group.members)
        for ci, host, own in branches_at.get(identity, ()):
            text += branch_text(ci, host, own)
        return text + "".join(suffixes[identity])

    texts, segments = [], []
    for ci, (chain, _) in enumerate(chains):
        if ci in attachments:
            continue
        text = "-".join(token(identity) for identity in chain)
        if ci in cyclic_tags:
            tag = cyclic_tags[ci]
            text = f"{tag}-{text}-{tag}"
            generated_markers.extend(
                (
                    LayoutMarker(Endpoint(chain[0], 1), tag, "terminal"),
                    LayoutMarker(Endpoint(chain[-1], 2), tag, "terminal"),
                )
            )
        texts.append(text)
        members = set(chain)
        for group in groups:
            if group.host in members:
                members.update(group.members)
        segments.append(
            LayoutSegment(len(segments), chain, tuple(i for i in nodes if i in members))
        )
    order = tuple(explicit_order + pendant_order)
    if len(order) != len(nodes) or set(order) != set(nodes):
        raise ValueError("Serialization did not emit each occurrence exactly once")
    preserved_markers = tuple(
        replace(marker, span=None)
        for marker in (peptide.layout.markers if peptide.layout else ())
        if marker.endpoint in represented
    )
    text = "%".join(texts)
    layout = NotationLayout(
        None,
        tuple(segments),
        tuple(groups),
        preserved_markers + tuple(generated_markers),
    )
    return Serialization(text, order, layout)
