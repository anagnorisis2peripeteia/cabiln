"""Library-free notation rules shared by parsing, editing, and serialization.

Legacy text adapters decide which inputs to transform. This module emits a chosen
chain layout and preserves tracked source when normalizing legacy bracket arms;
it neither resolves monomers nor chooses a chemical interpretation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from pyPept.source import group as _source_group
from pyPept.source import join as _source_join

_ENTRY = r"((?:[A-Za-z]\w*|!\w+))\((\d+),(\d+)\)"
BRACKET_ENTRY_RE = re.compile(r"\." + _ENTRY)
MIXED_ENTRY_RE = re.compile(r"\[\.([^\[\]]*)\]|\." + _ENTRY)
_CHAIN_TOKEN_RE = re.compile(r"([A-Za-z_]\w*)(?:\((\d+),(\d+)\))?")
_CROSSLINK = r"(!\w+)(?:\((\d+),(\d+)\))?"
INLINE_BOND_RE = re.compile(r"\." + _CROSSLINK)
_TEXT_CROSSLINK_RE = re.compile(r"<[^>]*>|" + _CROSSLINK)


@dataclass(frozen=True)
class BracketEntry:
    """One fully consumed entry, retaining its original source slices."""

    text: str
    token: str
    previous_slot: str
    own_slot: str


@dataclass(frozen=True)
class BracketArm:
    text: str
    entries: tuple[BracketEntry, ...]


@dataclass(frozen=True)
class ChainEntry:
    token: str
    slots: tuple[str, str] | None
    markers: tuple[tuple[str, tuple[str, str] | None], ...]


def parse_chain_entry(text):
    """Read a complete legacy branch token; return None for unsupported text.

    The optional pair on the monomer describes a chain continuation. Each
    crosslink has its own optional pair. Keeping those annotations separate
    prevents a text formatter from discarding unknown parenthesized content.
    """
    match = _CHAIN_TOKEN_RE.match(text)
    if match is None:
        return None
    token = _source_group(match, 1)
    slots = (match.group(2), match.group(3)) if match.group(2) else None
    markers = []
    position = match.end()
    while position < len(text):
        marker = INLINE_BOND_RE.match(text, position)
        if marker is None:
            return None
        pair = (marker.group(2), marker.group(3)) if marker.group(2) else None
        markers.append((marker.group(1), pair))
        position = marker.end()
    return ChainEntry(token, slots, tuple(markers))


def parse_bracket_entries(content, *, allow_arms=True, allow_empty=False):
    """Read dot-prefixed entries without skipping prefixes, gaps or suffixes.

    A modern arm is a flat sequence branching from the current fragment. Legacy
    lowering also uses this grammar before it moves any entries into arms.
    The slices retain SourceText origins; parsing itself records no occurrences.
    """
    entries = []
    position = 0
    pattern = MIXED_ENTRY_RE if allow_arms else BRACKET_ENTRY_RE
    while position < len(content):
        match = pattern.match(content, position)
        if match is None:
            reason = (
                "has unrecognised content" if entries else "contains no valid entries"
            )
            raise ValueError(
                f"Sequential {'sub-bracket' if not allow_arms else 'bracket'} "
                f"{reason} at offset {position}: "
                f"{content[position:position + 40]!r}; expected .Entry(r,r)"
                + (" or [.arm(r,r)]." if allow_arms else ".")
            )
        text = _source_group(match)
        if allow_arms and match.group(1) is not None:
            if not entries:
                raise ValueError("Sequential bracket must start with a flat entry.")
            arm = parse_bracket_entries("." + _source_group(match, 1), allow_arms=False)
            entries.append(BracketArm(text, arm))
        else:
            token_group = 2 if allow_arms else 1
            entries.append(BracketEntry(
                text, _source_group(match, token_group),
                match.group(token_group + 1), match.group(token_group + 2),
            ))
        position = match.end()
    if not entries and not allow_empty:
        raise ValueError("Sequential bracket has unrecognised content: no entries.")
    return tuple(entries)


def _bracket_end(text, start):
    """Return the end of one balanced bracket, rejecting mixed delimiters."""
    text = str(text)
    stack = []
    for index in range(start, len(text)):
        character = text[index]
        if character in "[{":
            stack.append("]" if character == "[" else "}")
        elif character in "]}":
            if not stack or character != stack.pop():
                raise ValueError(
                    f"Sequential bracket has mismatched delimiter at offset {index}."
                )
            if not stack:
                return index + 1
    raise ValueError(f"Sequential bracket at offset {start} is not closed.")


def bracket_regions(text):
    """Yield complete .[...] / .{...} spans, including any nested arms."""
    text = str(text)
    position = 0
    while position < len(text) - 1:
        if text[position] == "." and text[position + 1] in "[{":
            end = _bracket_end(text, position + 1)
            yield position, end
            position = end
        else:
            position += 1


def crosslink_declarations(entries):
    """Yield declared slot pairs in bracket traversal order, including arms."""
    for entry in entries:
        if isinstance(entry, BracketArm):
            yield from crosslink_declarations(entry.entries)
        elif entry.token.startswith("!"):
            yield entry.token, entry.previous_slot, entry.own_slot


def text_crosslink_declarations(text):
    """Read explicit marker facts for legacy text-conversion guards.

    A marker has the same declaration inside or outside a bracket. This reader
    does not claim the surrounding text is valid; adapters still consume their
    supported branch grammar completely. Synthetic SMILES remain opaque.
    """
    for marker in _TEXT_CROSSLINK_RE.finditer(text):
        if marker.group(2) is not None:
            yield marker.group(1), marker.group(2), marker.group(3)


def validate_crosslink_declarations(declarations, *, strict=True):
    """Resolve explicit endpoint pairs once, independent of their spelling.

    Parsing raises on conflicts. Legacy text adapters use strict=False and
    preserve the input if this returns None, without attempting a conversion.
    """
    seen = {}
    for token, previous, own in declarations:
        if token not in seen:
            seen[token] = (previous, own)
        else:
            first_previous, first_own = seen[token]
            if (previous, own) != (first_own, first_previous):
                if not strict:
                    return None
                raise ValueError(
                    f"Bond {token!r}: R-group conflict. "
                    f"First endpoint declared R{first_previous}→partner R{first_own}; "
                    f"second endpoint declares R{previous}→partner R{own} "
                    f"(expected inverse R{first_own}→partner R{first_previous}).")
    return seen


def supports_bracket_token(token):
    """Whether the accepted bracket grammar can spell this exact token."""
    return BRACKET_ENTRY_RE.fullmatch(f".{token}(1,1)") is not None


def legacy_attachment_slot(slot: int) -> int:
    """Legacy BILN/HELM R3 names the sidechain, stored at CABILN R4."""
    return 4 if slot == 3 else slot


@dataclass(frozen=True)
class BracketEmission:
    text: str
    order: tuple[int, ...]
    arm: tuple[int, ...]


def bracket_chain(entries, anchor, host_slot, own_slot, suffixes=None):
    """Write a chain attached at its end or midpoint, without resolving symbols.

    Entries are (token, incoming slot, outgoing slot) in the explicit chain's
    order. Missing slots mean R2-to-R1. Entries before the attachment read in
    reverse; a midpoint attachment gives the forward continuation its own arm.
    Returned indices record the exact order chosen here for both producers.
    """
    after = tuple(range(anchor + 1, len(entries)))
    before = tuple(range(anchor - 1, -1, -1))
    order = (anchor,) + after + before
    suffixes = suffixes or {}

    def entry(index, slots):
        return f"{entries[index][0]}({slots[0]},{slots[1]})" + suffixes.get(index, "")

    forward = []
    for index in after:
        _, previous, own = entries[index]
        forward.append(entry(index, (previous, own) if previous and own else (2, 1)))
    backward = []
    for index in before:
        _, previous, own = entries[index]
        backward.append(entry(index, (own, previous) if previous and own else (1, 2)))
    text = entry(anchor, (host_slot, own_slot))
    if forward:
        continuation = ".".join(forward)
        text += "[." + continuation + "]" if backward else "." + continuation
    if backward:
        text += "." + ".".join(backward)
    return BracketEmission(".[" + text + "]", order, after if after and before else ())


def normalize_legacy_brackets(seg):
    """Convert old ``[[A.B].C]`` nested notation to new ``[A[.B][.C]]`` form.

    Old notation used nesting to express multi-arm hubs: ``[[TBMB.C].!3]``
    meant both C and !3 bond to TBMB.  New notation uses explicit sub-bracket
    arms: ``[TBMB[.C][.!3]]``.  This function converts the old form on the way
    into the parser so both notations assemble identically.
    """
    out = []
    previous_end = 0
    for start, end in bracket_regions(seg):
        if seg.startswith(".[[", start):
            out.append(seg[previous_end : start + 1])
            out.append(_flatten_one_nested(seg[start + 1 : end]))
            previous_end = end
    if not out:
        return seg
    out.append(seg[previous_end:])
    return _source_join("", out)


def _flatten_one_nested(s):
    """Convert a single nested bracket to new sub-bracket arm notation.

    Input ``s`` is the bracket expression starting at the outer ``[``.
    Returns ``[anchor[.inner_rest...][.outer...]]`` — every arm from the
    anchor is an explicit ``[.Entry(r,r)]`` sub-bracket so the parser's
    pointer semantics are unambiguous.
    """
    if not s.startswith("[[") or _bracket_end(s, 0) != len(s):
        raise ValueError(
            "Legacy nested bracket must be one complete [[...].Entry(r,r)] group."
        )
    inner_end = _bracket_end(s, 1)
    inner_content = s[2 : inner_end - 1]
    outer_part = s[inner_end:-1]
    inner_entries = parse_bracket_entries("." + inner_content)
    outer_entries = parse_bracket_entries(
        outer_part, allow_arms=False, allow_empty=True
    )

    def arms(entries):
        # Drop only the leading dot. Tokens and slot declarations keep their
        # original slices, including ownership spans in tracked input.
        return _source_join("", ("[." + entry.text[1:] + "]" for entry in entries))

    # If inner content already uses sub-bracket notation, preserve it and just
    # append outer entries as additional arms.
    if any(isinstance(entry, BracketArm) for entry in inner_entries):
        return s[:1] + inner_content + arms(outer_entries) + s[-1:]

    # Legacy flat inner content: split anchor from remaining chain entries.
    anchor = inner_entries[0].text[1:]
    return s[:1] + anchor + arms(inner_entries[1:]) + arms(outer_entries) + s[-1:]
