"""Source-preserving edits verified against resolved peptide definitions.

A selection belongs to one source revision. An edit changes only the selected
source occurrences, reparses them, and verifies the requested change in slot
connections before assembling the product. No second chemical graph is kept.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256

from pyPept.molecule import Molecule
from pyPept.attachments import reaction_for_types
from pyPept.notation import normalize_legacy_brackets, supports_bracket_token
from pyPept.peptide import Connection, Endpoint, Peptide, serialize
from pyPept.sequence import Sequence
from pyPept.source import SourceText, Span, join, origin_span


@dataclass(frozen=True)
class Selection:
    revision: str
    index: int


@dataclass(frozen=True)
class _Edit:
    span: Span
    replacement: str


class EditError(ValueError):
    """An edit cannot preserve the requested source or assembly topology."""


class PeptideDocument:
    """One source revision and its derived, currently configured Sequence."""

    def __init__(self, source: str):
        self.source = source
        self.revision = sha256(source.encode()).hexdigest()
        self.sequence = Sequence(source, track_source=True)
        self.peptide = Peptide.from_sequence(self.sequence)
        self._text = SourceText.original(source)

    def select(self, index: int) -> Selection:
        if not 0 <= index < len(self.sequence.s_monomers):
            raise ValueError(f"Residue {index} does not exist")
        return Selection(self.revision, index)

    def _index(self, selection: Selection) -> int:
        if selection.revision != self.revision:
            raise EditError(
                "Selection belongs to an older source revision; select again"
            )
        return self.select(selection.index).index

    def _editable_brackets(self, *indices):
        # Legacy [[anchor.arm].arm] lowers sibling entries into generated arms.
        # Extending such an entry in the old spelling would attach to its hub.
        # Use that SAME lowering helper on the selected bracket first, checking
        # that node order, templates and connections remain identical.
        brackets = set()
        for index in indices:
            occurrence = self.sequence.s_sources[index]
            if occurrence.legacy_arm:
                brackets.add(occurrence.bracket)
        if not brackets:
            return self
        source = self.source
        for span in sorted(brackets, reverse=True):
            source = (
                source[: span.start]
                + normalize_legacy_brackets(source[span.start : span.end])
                + source[span.end :]
            )
        updated = PeptideDocument(source)
        if (
            updated.sequence.s_biln != self.sequence.s_biln
            or updated.peptide.connections != self.peptide.connections
            or [node.definition for node in updated.peptide.occurrences]
            != [node.definition for node in self.peptide.occurrences]
        ):
            raise EditError("Legacy bracket normalization would change the sequence")
        return updated

    def _free_slot(self, index: int, slot: int):
        endpoint = Endpoint(index, slot)
        try:
            self.peptide.site(endpoint)
        except ValueError as error:
            raise ValueError(f"Residue {index} has no R{slot} attachment") from error
        if self.peptide.connection_at(endpoint) is not None:
            raise ValueError(f"Residue {index} R{slot} is already bonded")

    def replacement_requirements(self, selection: Selection):
        """Describe every occupied site, including branches and ring closures."""
        index = self._index(selection)
        requirements = []
        for edge in self.peptide.connections:
            for own, other in (edge.endpoints, tuple(reversed(edge.endpoints))):
                if own.occurrence_id != index:
                    continue
                requirements.append({
                    "slot": own.slot,
                    "chem_type": self.peptide.site(own).chem_type,
                    "partner_idx": other.occurrence_id,
                    "partner_abbr": self.peptide.occurrence(other.occurrence_id).symbol,
                    "partner_slot": other.slot,
                    "partner_chem_type": self.peptide.site(other).chem_type,
                    "internal": other.occurrence_id == index,
                })
        return sorted(requirements, key=lambda item: item["slot"])

    @staticmethod
    def replacement_options(requirements, sites):
        """Find distinct compatible sites for every existing connection.

        Options describe detected chemistry, not a guarantee that a reaction's
        structural pattern will assemble. Applying the chosen mapping validates
        the complete product. Equal slot numbers and chemistry are preferred;
        the caller still presents the mapping for review.
        """
        types = {site["slot"]: site["chem_type"] for site in sites}
        choices = {
            need["slot"]: [site["slot"] for site in sorted(
                sites,
                key=lambda site: (
                    site["slot"] != need["slot"],
                    site["chem_type"] != need["chem_type"], site["slot"],
                ),
            ) if (
                any(other != site["slot"] and reaction_for_types(site["chem_type"], chemistry)
                    for other, chemistry in types.items())
                if need.get("internal")
                else reaction_for_types(site["chem_type"], need["partner_chem_type"])
            )]
            for need in requirements
        }

        def match(fixed):
            owners = {target: slot for slot, target in fixed.items()}

            def assign(slot, visited):
                for target in choices[slot]:
                    if target in visited or owners.get(target) in fixed:
                        continue
                    visited.add(target)
                    if target not in owners or assign(owners[target], visited):
                        owners[target] = slot
                        return True
                return False

            remaining = sorted(set(choices) - set(fixed), key=lambda s: len(choices[s]))
            if all(assign(slot, set()) for slot in remaining):
                return {slot: target for target, slot in owners.items()}
            return None

        internal = [(need["slot"], need["partner_slot"]) for need in requirements
                    if need.get("internal") and need["slot"] < need["partner_slot"]]

        def match_internal(pairs, fixed):
            # A closure wholly inside the replaced monomer depends on BOTH new
            # sites. Fix compatible pairs before matching its external bonds.
            if not pairs:
                return match(fixed)
            left, right = pairs[0]
            for a in choices[left]:
                for b in choices[right]:
                    if (a == b or a in fixed.values() or b in fixed.values()
                            or not reaction_for_types(types[a], types[b])):
                        continue
                    result = match_internal(pairs[1:], {**fixed, left: a, right: b})
                    if result is not None:
                        return result
            return None

        mapping = match_internal(internal, {})
        if mapping is None:
            return None
        return {"choices": choices, "mapping": mapping}

    def replace_monomer(self, selection: Selection, symbol: str, mapping: dict[int, int]) -> str:
        """Replace one definition, preserving all neighbours and explicit site choices."""
        index = self._index(selection)
        replacement = Peptide.from_sequence(Sequence(symbol))
        if len(replacement.occurrences) != 1 or replacement.connections:
            raise EditError("Choose a single replacement monomer")
        node = replacement.occurrences[0]
        sites = [{"slot": s.slot, "chem_type": s.chem_type} for s in node.sites]
        requirements = self.replacement_requirements(selection)
        options = self.replacement_options(requirements, sites)
        if options is None:
            raise EditError("This monomer cannot preserve all existing connections")
        if set(mapping) != set(options["choices"]):
            raise EditError("Map every occupied R-group, and only occupied R-groups")
        if len(set(mapping.values())) != len(mapping):
            raise EditError("Each existing connection needs a different replacement R-group")
        for old, new in mapping.items():
            if new not in options["choices"][old]:
                raise EditError(f"Replacement R{new} is incompatible with the neighbour at R{old}")
        types = {site.slot: site.chem_type for site in node.sites}
        for need in requirements:
            if need["internal"] and not reaction_for_types(
                types[mapping[need["slot"]]], types[mapping[need["partner_slot"]]]
            ):
                raise EditError("The mapped sites cannot close the replacement monomer's internal connection")

        def endpoint(site):
            return Endpoint(index, mapping[site.slot]) if site.occurrence_id == index else site

        nodes = tuple(
            replace(node, id=index, token=node.definition.token)
            if old.id == index else old
            for old in self.peptide.occurrences
        )
        expected = Peptide(nodes, tuple(
            Connection(endpoint(edge.left), endpoint(edge.right), edge.label)
            for edge in self.peptide.connections
        ))

        # The common same-number swap preserves spelling, brackets and spacing.
        # Renumbered sites are written from the graph, never by guessing which
        # parenthesised numbers in nested notation refer to this occurrence.
        text = None
        if all(old == new for old, new in mapping.items()):
            span = self.sequence.s_sources[index].token
            candidate = self.source[:span.start] + symbol + self.source[span.end:]
            try:
                parsed = Peptide.from_sequence(Sequence(candidate, track_source=True))
                if self._same_replacement(parsed, expected):
                    text = candidate
            except ValueError:
                pass  # A library symbol can require an explicit notation segment.
        if text is None:
            notation = "bracket" if self.peptide.layout.groups else "percent"
            emitted = serialize(replace(expected, occurrences=tuple(
                replace(item, token=item.definition.token) for item in expected.occurrences
            )), notation)
            text = emitted.text
            parsed = Peptide.from_sequence(Sequence(text, track_source=True), emitted.occurrence_order)
            if not self._same_replacement(parsed, expected):
                raise EditError("Replacement notation would change another monomer or connection")
        # Validate both the requested graph and its emitted representation.
        from pyPept.structure import compare_structures

        intended = Molecule(expected, depiction=None).get_molecule("ROMol")
        actual = Molecule(parsed, depiction=None).get_molecule("ROMol")
        if not compare_structures(intended, actual).exact:
            raise EditError("Replacement notation would change the assembled product")
        return text

    @staticmethod
    def _same_replacement(actual, expected):
        return (
            {node.id: node.definition for node in actual.occurrences}
            == {node.id: node.definition for node in expected.occurrences}
            and set(actual.connections) == set(expected.connections)
        )

    @staticmethod
    def _new_definition(symbol: str, *slots: int):
        sequence = Sequence(symbol)
        if len(sequence.s_monomers) != 1:
            raise ValueError("Choose a single monomer to insert")
        peptide = Peptide.from_sequence(sequence)
        for slot in slots:
            try:
                peptide.site(Endpoint(0, slot))
            except ValueError as error:
                raise ValueError(
                    f"Monomer {symbol!r} has no R{slot} attachment"
                ) from error
        return peptide.occurrences[0].definition

    def _tag(self):
        used = {edge.label for edge in self.peptide.connections}
        number = 1
        while f"!{number}" in used:
            number += 1
        return f"!{number}"

    def _append_to_occurrence(self, index: int, suffix: str) -> list[_Edit]:
        occurrence = self.sequence.s_sources[index]
        if occurrence.kind == "inline":
            # Give this attachment its own continuation without changing the
            # backbone host's pointer. Copy its original token and slot text.
            body = self._text[occurrence.entry.start + 1 : occurrence.entry.end]
            if occurrence.bracketable:
                return [_Edit(occurrence.entry, ".[" + body + suffix + "]")]
            # Some library symbols (for example leading underscores) are legal
            # inline or explicit tokens but cannot occur in a bracket. Move
            # only this occurrence into an explicit segment, preserving its
            # original characters and existing attachment.
            bond = next(
                edge
                for edge in self.peptide.connections
                if any(endpoint.occurrence_id == index for endpoint in edge.endpoints)
            )
            own = next(
                endpoint
                for endpoint in bond.endpoints
                if endpoint.occurrence_id == index
            )
            host = next(endpoint for endpoint in bond.endpoints if endpoint != own)
            own_slot, host_slot = own.slot, host.slot
            tag = f"!source{index}"
            while tag in {edge.label for edge in self.peptide.connections}:
                tag += "_"
            token = self._text[occurrence.token.start : occurrence.token.end]
            insertion = len(self.source.rstrip())
            return [
                _Edit(occurrence.entry, f".{tag}({host_slot},{own_slot})"),
                _Edit(
                    Span(insertion, insertion),
                    "%" + token + f".{tag}({own_slot},{host_slot})" + suffix,
                ),
            ]
        return [_Edit(Span(occurrence.entry.end, occurrence.entry.end), suffix)]

    def connect(
        self, host: Selection, host_slot: int, target: Selection, target_slot: int
    ) -> str:
        a, b = self._index(host), self._index(target)
        editable = self._editable_brackets(a, b)
        if editable is not self:
            return editable.connect(
                editable.select(a), host_slot, editable.select(b), target_slot
            )
        if a == b:
            raise ValueError("Choose two different residues for a crosslink")
        self._free_slot(a, host_slot)
        self._free_slot(b, target_slot)
        tag = self._tag()
        edits = self._append_to_occurrence(
            a, f".{tag}({host_slot},{target_slot})"
        ) + self._append_to_occurrence(b, f".{tag}({target_slot},{host_slot})")
        expected = set(self.peptide.connections) | {
            Connection(Endpoint(a, host_slot), Endpoint(b, target_slot))
        }
        return self._apply(edits, expected)

    def attach(
        self, host: Selection, host_slot: int, symbol: str, new_slot: int
    ) -> str:
        index = self._index(host)
        editable = self._editable_brackets(index)
        if editable is not self:
            return editable.attach(editable.select(index), host_slot, symbol, new_slot)
        self._free_slot(index, host_slot)
        definition = self._new_definition(symbol, new_slot)
        occurrence = self.sequence.s_sources[index]
        chain = self.sequence.s_chains["s_monomerIDs"][
            self.sequence.s_monomers[index]["m_chainID"]
        ]
        start, end = occurrence.entry.start, occurrence.entry.end

        bracket_token = supports_bracket_token(symbol)
        if occurrence.kind == "explicit" and (
            bracket_token or (host_slot, new_slot) in ((2, 1), (1, 2))
        ):
            if (host_slot, new_slot) == (2, 1) and index == chain[-1]:
                edit = _Edit(Span(end, end), "-" + symbol)
            elif (host_slot, new_slot) == (1, 2) and index == chain[0]:
                edit = _Edit(Span(start, start), symbol + "-")
            else:
                edit = _Edit(Span(end, end), f".{{{symbol}({host_slot},{new_slot})}}")
            edits = [edit]
        elif bracket_token and (occurrence.kind == "inline" or occurrence.terminal):
            edits = self._append_to_occurrence(
                index, f".{symbol}({host_slot},{new_slot})"
            )
        elif bracket_token:
            # A child scope returns to this occurrence before any existing
            # children or flat continuation, so their attachment host survives.
            edits = self._append_to_occurrence(
                index, f".{{{symbol}({host_slot},{new_slot})}}"
            )
        else:
            # A crosslink annotation leaves the parser's current pointer intact,
            # preserving any following flat entries or nested arms. Put the new
            # monomer in an explicit segment instead of redirecting those edges.
            tag = self._tag()
            insertion = len(self.source.rstrip())
            edits = self._append_to_occurrence(index, f".{tag}({host_slot},{new_slot})")
            edits.append(
                _Edit(
                    Span(insertion, insertion),
                    f"%{symbol}.{tag}({new_slot},{host_slot})",
                )
            )
        new_index = len(self.sequence.s_monomers)
        expected = set(self.peptide.connections) | {
            Connection(Endpoint(index, host_slot), Endpoint(new_index, new_slot))
        }
        return self._apply(edits, expected, definition)

    def insert_backbone(self, after: Selection, symbol: str) -> str:
        index = self._index(after)
        occurrence = self.sequence.s_sources[index]
        if occurrence.kind != "explicit":
            raise EditError(
                "Backbone insertion requires a residue in an explicit chain"
            )
        expected = set(self.peptide.connections)
        endpoint = Endpoint(index, 2)
        try:
            outgoing = self.peptide.connection_at(endpoint)
        except ValueError as error:
            raise ValueError(f"Residue {index} has no R2 attachment") from error
        definition = self._new_definition(symbol, 1, *([2] if outgoing else []))
        new_index = len(self.sequence.s_monomers)
        if outgoing is not None:
            other = next(site for site in outgoing.endpoints if site != endpoint)
            expected.remove(outgoing)
            expected.add(Connection(Endpoint(new_index, 2), other))
        expected.add(Connection(endpoint, Endpoint(new_index, 1)))
        end = occurrence.entry.end
        suffix = "-" + symbol
        edits = []
        if outgoing is not None:
            marker = next(
                (
                    m
                    for m in self.peptide.layout.markers
                    if m.endpoint == endpoint and m.kind != "terminal"
                ),
                None,
            )
            if marker is not None:
                removal, annotation = self._move_marker(marker)
                edits.append(removal)
                suffix += annotation
        edits.append(_Edit(Span(end, end), suffix))
        return self._apply(edits, expected, definition)

    def _move_marker(self, marker):
        """Extract one endpoint, retaining the spelling of any surrounding bracket."""
        if marker.group is None:
            return (
                _Edit(marker.span, ""),
                self._text[marker.span.start : marker.span.end],
            )
        groups = {group.id: group for group in self.peptide.layout.groups}
        group = groups[marker.group]

        def contains_only_marker(candidate):
            return not candidate.members and all(
                other == marker
                for other in self.peptide.layout.markers
                if candidate.span.start <= other.span.start
                and other.span.end <= candidate.span.end
            )

        # Remove scopes made empty by the move, including marker-only ancestors.
        # A group containing monomers or another endpoint must stay on its host.
        empty = None
        candidate = group
        while candidate is not None and contains_only_marker(candidate):
            empty = candidate
            candidate = groups.get(candidate.parent)
        if empty is not None:
            body = self._text[empty.span.start : empty.span.end]
            if empty.protected and empty.opening == "[":
                # Protection inherited from a parent must travel with the
                # endpoint when this child moves outside that parent.
                offset = 2 if body.startswith(".[") else 1
                if body.startswith((".[", "[")):
                    body = ".{" + body[offset:-1] + "}"
            if body.startswith((".[", ".{")):
                return _Edit(empty.span, ""), body
            if body.startswith(("[", "{")):
                return _Edit(empty.span, ""), "." + body

        # Only remove this marker, not everything before the next marker: that
        # interval can also contain monomers or independently hosted children.
        body = self._text[marker.span.start : marker.span.end]
        removed, replacement = marker.span, ""
        prefix = 2 if self.source[group.span.start] == "." else 1
        if marker.span.start == group.span.start + prefix:
            following = self.source[marker.span.end : marker.span.end + 2]
            if following.startswith(".") and following[1:] not in ("[", "{"):
                if not body.startswith("."):
                    removed = Span(marker.span.start, marker.span.end + 1)
            elif following.startswith(("[", "{")):
                # A child newly at the start must not become a legacy hub wrapper.
                replacement = "."
        if body.startswith("."):
            body = body[1:]
        opening, closing = ("{", "}") if group.protected else ("[", "]")
        annotation = "." + opening + body + closing
        return _Edit(removed, replacement), annotation

    def _apply(self, edits, expected_edges, new_definition=None):
        pieces, cursor = [], 0
        for edit in sorted(edits, key=lambda item: (item.span.start, item.span.end)):
            if edit.span.start < cursor:
                raise EditError("The selected source edits overlap")
            pieces.extend((self._text[cursor : edit.span.start], edit.replacement))
            cursor = edit.span.end
        pieces.append(self._text[cursor:])
        mapped_result = join("", pieces)
        result = str(mapped_result)
        updated = Sequence(result, track_source=True)

        # Original character origins identify each retained instance even if the
        # parser moves appended branches or the same symbol occurs many times.
        original_ids = {s.token: i for i, s in enumerate(self.sequence.s_sources)}
        identities, seen = {}, set()
        for index, occurrence in enumerate(updated.s_sources):
            old_span = origin_span(
                mapped_result[occurrence.token.start : occurrence.token.end]
            )
            if old_span is None:
                if new_definition is None or len(self.sequence.s_monomers) in seen:
                    raise EditError("Edit created an unexpected monomer")
                identity = len(self.sequence.s_monomers)
            else:
                identity = original_ids.get(old_span)
                if identity is None or identity in seen:
                    raise EditError(
                        "Edit changed or duplicated an existing source occurrence"
                    )
            identities[index] = identity
            seen.add(identity)
        definitions = [node.definition for node in self.peptide.occurrences]
        if new_definition is not None:
            definitions.append(new_definition)
        expected_count = len(definitions)
        if seen != set(range(expected_count)):
            raise EditError("Edit removed an existing monomer")
        updated_peptide = Peptide.from_sequence(
            updated, tuple(identities[i] for i in range(len(identities)))
        )
        if any(
            node.definition != definitions[node.id]
            for node in updated_peptide.occurrences
        ):
            raise EditError("Edit changed the identity or attachments of a monomer")
        if set(updated_peptide.connections) != expected_edges:
            raise EditError(
                "Edit would change connections beyond the selected attachment"
            )
        Molecule(updated_peptide)
        return result
