"""Source-preserving edits over the existing Sequence assembly model.

A selection belongs to one source revision. An edit changes only the selected
source occurrences, reparses them, and verifies the requested change in slot
connections before assembling the product. No second chemical graph is kept.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256

from rdkit import Chem

from pyPept.molecule import Molecule
from pyPept.sequence import (
    _BRACKET_ENTRY_RE,
    Sequence,
    _flatten_nested_brackets,
    _rgroup_atom_idx,
)
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


def _edge(a, slot_a, b, slot_b):
    return tuple(sorted(((a, slot_a), (b, slot_b))))


def _edges(sequence, identities=None):
    if identities is None:
        identities = {i: i for i in range(len(sequence.s_monomers))}
    return {
        _edge(identities[b[0]], b[4], identities[b[2]], b[5]) for b in sequence.s_bonds
    }


def _template(monomer):
    return (
        monomer["m_abbr"],
        Chem.MolToSmiles(monomer["m_romol"]),
        tuple(monomer["m_Rgroups"]),
    )


class PeptideDocument:
    """One source revision and its derived, currently configured Sequence."""

    def __init__(self, source: str):
        self.source = source
        self.revision = sha256(source.encode()).hexdigest()
        self.sequence = Sequence(source, track_source=True)
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
            if occurrence.arm is not None and not self.source[
                occurrence.arm.start : occurrence.arm.end
            ].startswith("["):
                brackets.add(occurrence.bracket)
        if not brackets:
            return self
        source = self.source
        for span in sorted(brackets, reverse=True):
            source = (
                source[: span.start]
                + _flatten_nested_brackets(source[span.start : span.end])
                + source[span.end :]
            )
        updated = PeptideDocument(source)
        if (
            updated.sequence.s_biln != self.sequence.s_biln
            or updated.sequence.s_bonds != self.sequence.s_bonds
            or [_template(m) for m in updated.sequence.s_monomers]
            != [_template(m) for m in self.sequence.s_monomers]
        ):
            raise EditError("Legacy bracket normalization would change the sequence")
        return updated

    def _free_slot(self, index: int, slot: int):
        monomer = self.sequence.s_monomers[index]
        if _rgroup_atom_idx(monomer["m_romol"], slot) is None:
            raise ValueError(f"Residue {index} has no R{slot} attachment")
        if any((index, slot) in edge for edge in _edges(self.sequence)):
            raise ValueError(f"Residue {index} R{slot} is already bonded")

    @staticmethod
    def _new_monomer(symbol: str, *slots: int):
        sequence = Sequence(symbol)
        if len(sequence.s_monomers) != 1:
            raise ValueError("Choose a single monomer to insert")
        monomer = sequence.s_monomers[0]
        for slot in slots:
            if _rgroup_atom_idx(monomer["m_romol"], slot) is None:
                raise ValueError(f"Monomer {symbol!r} has no R{slot} attachment")
        return monomer

    def _tag(self):
        used = self.sequence.s_biln.tracker.labels
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
            if _BRACKET_ENTRY_RE.fullmatch("." + body):
                return [_Edit(occurrence.entry, ".[" + body + suffix + "]")]
            # Some library symbols (for example leading underscores) are legal
            # inline or explicit tokens but cannot occur in a bracket. Move
            # only this occurrence into an explicit segment, preserving its
            # original characters and existing attachment.
            bond = next(b for b in self.sequence.s_bonds if index in (b[0], b[2]))
            own_slot, host_slot = (
                (bond[4], bond[5]) if bond[0] == index else (bond[5], bond[4])
            )
            tag = f"!source{index}"
            while tag in self.sequence.s_biln.tracker.labels:
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
        expected = _edges(self.sequence) | {_edge(a, host_slot, b, target_slot)}
        return self._apply(edits, expected)

    def attach(
        self, host: Selection, host_slot: int, symbol: str, new_slot: int
    ) -> str:
        index = self._index(host)
        editable = self._editable_brackets(index)
        if editable is not self:
            return editable.attach(editable.select(index), host_slot, symbol, new_slot)
        self._free_slot(index, host_slot)
        monomer = self._new_monomer(symbol, new_slot)
        occurrence = self.sequence.s_sources[index]
        chain = self.sequence.s_chains["s_monomerIDs"][
            self.sequence.s_monomers[index]["m_chainID"]
        ]
        start, end = occurrence.entry.start, occurrence.entry.end

        bracket_token = bool(
            _BRACKET_ENTRY_RE.fullmatch(f".{symbol}({host_slot},{new_slot})")
        )
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
        expected = _edges(self.sequence) | {
            _edge(index, host_slot, new_index, new_slot)
        }
        return self._apply(edits, expected, monomer)

    def insert_backbone(self, after: Selection, symbol: str) -> str:
        index = self._index(after)
        occurrence = self.sequence.s_sources[index]
        if occurrence.kind != "explicit":
            raise EditError(
                "Backbone insertion requires a residue in an explicit chain"
            )
        expected = _edges(self.sequence)
        outgoing = next((edge for edge in expected if (index, 2) in edge), None)
        monomer = self._new_monomer(symbol, 1, *([2] if outgoing else []))
        if _rgroup_atom_idx(self.sequence.s_monomers[index]["m_romol"], 2) is None:
            raise ValueError(f"Residue {index} has no R2 attachment")
        new_index = len(self.sequence.s_monomers)
        if outgoing is not None:
            other, other_slot = next(
                endpoint for endpoint in outgoing if endpoint != (index, 2)
            )
            expected.remove(outgoing)
            expected.add(_edge(new_index, 2, other, other_slot))
        expected.add(_edge(index, 2, new_index, 1))
        end = occurrence.entry.end
        suffix = "-" + symbol
        edits = []
        if outgoing is not None:
            marker = next(
                (
                    m
                    for m in self.sequence.s_biln.tracker.markers
                    if m.slot == 2
                    and m.owner == occurrence.entry
                    and m.kind != "terminal"
                ),
                None,
            )
            if marker is not None:
                removed, annotation = self._move_marker(marker)
                edits.append(_Edit(removed, ""))
                suffix += annotation
        edits.append(_Edit(Span(end, end), suffix))
        return self._apply(edits, expected, monomer)

    def _move_marker(self, marker):
        """Extract one endpoint, retaining the spelling of any surrounding bracket."""
        if marker.bracket is None:
            return marker.span, self._text[marker.span.start : marker.span.end]
        siblings = sorted(
            (
                m
                for m in self.sequence.s_biln.tracker.markers
                if m.bracket == marker.bracket
            ),
            key=lambda m: m.span.start,
        )
        bracket = marker.bracket
        if len(siblings) == 1:
            return bracket, self._text[bracket.start : bracket.end]
        body = self._text[marker.span.start : marker.span.end]
        if marker == siblings[0]:
            # The first entry has no leading dot. Remove the separator before
            # the next recorded entry so its name becomes the first entry.
            removed = Span(marker.span.start, siblings[1].span.start + 1)
        else:
            removed = marker.span
            body = body[1:]
        annotation = (
            self._text[bracket.start : bracket.start + 2]
            + body
            + self._text[bracket.end - 1 : bracket.end]
        )
        return removed, annotation

    def _apply(self, edits, expected_edges, new_monomer=None):
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
                if new_monomer is None or len(self.sequence.s_monomers) in seen:
                    raise EditError("Edit created an unexpected monomer")
                identity = len(self.sequence.s_monomers)
                template = new_monomer
            else:
                identity = original_ids.get(old_span)
                if identity is None or identity in seen:
                    raise EditError(
                        "Edit changed or duplicated an existing source occurrence"
                    )
                template = self.sequence.s_monomers[identity]
            if _template(updated.s_monomers[index]) != _template(template):
                raise EditError("Edit changed the identity or attachments of a monomer")
            identities[index] = identity
            seen.add(identity)
        expected_count = len(self.sequence.s_monomers) + (new_monomer is not None)
        if seen != set(range(expected_count)):
            raise EditError("Edit removed an existing monomer")
        if _edges(updated, identities) != expected_edges:
            raise EditError(
                "Edit would change connections beyond the selected attachment"
            )
        Molecule(updated)
        return result
