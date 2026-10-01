"""Library-free notation rules shared by parsing, editing, and serialization.

Legacy text adapters decide which inputs to transform. This module emits a chosen
chain layout and preserves tracked source when normalizing legacy bracket arms;
it neither resolves monomers nor chooses a chemical interpretation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace

from pyPept.source import group as _source_group
from pyPept.source import join as _source_join

_MONOMER = (
    r"(?:_<[^<>]*>|<[^<>]*>_?|"
    r"[A-Za-z0-9_]\w*(?:\([A-Za-z_]\w*(?:,[A-Za-z_]\w*)*\))*)"
)
_ENTRY = r"(" + _MONOMER + r"|!\w+)\((\d+),(\d+)\)"
BRACKET_ENTRY_RE = re.compile(r"\." + _ENTRY)
_SCOPED_ENTRY_RE = re.compile(r"(" + _MONOMER + r"|!\w+)(?:\((\d+),(\d+)\))?")
_CHAIN_TOKEN_RE = re.compile(r"([A-Za-z_]\w*)(?:\((\d+),(\d+)\))?")
_CROSSLINK = r"(!\w+)(?:\((\d+),(\d+)\))?"
INLINE_BOND_RE = re.compile(r"\." + _CROSSLINK)
_TEXT_CROSSLINK_RE = re.compile(r"<[^>]*>|" + _CROSSLINK)

MAX_NOTATION_CHARACTERS = 1_000_000
MAX_BRACKET_DEPTH = 2048
MAX_BRACKET_ENTRIES = 100_000


@dataclass(frozen=True)
class BracketEntry:
    """One fully consumed entry, retaining its original source slices."""

    text: str
    token: str
    previous_slot: str | None
    own_slot: str | None


@dataclass(frozen=True)
class BracketArm:
    """One source scope; legacy arms may retain an entry span without brackets."""

    text: str
    entries: tuple[BracketEntry | BracketArm, ...]
    opening: str = "["
    closing: str = "]"
    legacy: bool = False
    normalized_legacy: bool = False
    # A retained protected legacy wrapper exports its hub, preserving outer arms.
    returns_anchor: bool = False


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
    """Read complete dot-prefixed entries through the recursive scope reader."""
    if not content and allow_empty:
        return ()
    if not content.startswith("."):
        raise ValueError("Sequential bracket must start with .Entry(r,r).")
    scope = parse_bracket_group("[" + content + "]")
    if not allow_arms and any(isinstance(e, BracketArm) for e in scope.entries):
        raise ValueError("Sequential sub-bracket does not allow another arm here.")
    return scope.entries


def _legacy_entries(entries):
    """Lower a leading legacy group into hub arms without rewriting its text."""
    wrapper = entries[0]
    inner = wrapper.entries
    if not inner or (isinstance(inner[0], BracketArm) and not inner[0].returns_anchor):
        raise ValueError("Legacy sequential bracket must have an anchor entry.")

    def arm(entry):
        return (
            entry
            if isinstance(entry, BracketArm)
            else BracketArm(entry.text, (entry,), legacy=True)
        )

    # The established spelling treats flat inner entries as independent arms.
    # An already explicit inner arm keeps the surrounding sequential semantics.
    if not any(isinstance(entry, BracketArm) for entry in inner):
        inner = (inner[0],) + tuple(arm(entry) for entry in inner[1:])
    if wrapper.opening == "{":
        # A protected legacy wrapper is meaningful source layout. Keep its
        # actual scope, and export its anchor to the outer legacy hub on close.
        inner = (
            replace(
                wrapper,
                entries=tuple(inner),
                returns_anchor=True,
                normalized_legacy=True,
            ),
        )
    return tuple(inner) + tuple(arm(entry) for entry in entries[1:])


def parse_bracket_group(text):
    """Read one complete scope with an explicit stack and original text slices.

    Modern ``.[...]``/``.{...}``, old ``[.arm]``, and legacy hub wrappers all
    produce the same entries. Every child has its own cursor in the consumer.
    No Python recursion or substring-rewrite loop determines parent ownership.
    """
    if len(text) > MAX_NOTATION_CHARACTERS:
        raise ValueError(
            f"Notation exceeds {MAX_NOTATION_CHARACTERS:,} characters; "
            "split it into smaller documents."
        )
    opening_at = 1 if text.startswith(".") else 0
    if opening_at >= len(text) or text[opening_at] not in "[{":
        raise ValueError("Sequential bracket must begin with '[' or '{'.")
    frames = [
        {
            "start": 0,
            "opening": str(text[opening_at]),
            "entries": [],
            "legacy": False,
        }
    ]
    position = opening_at + 1
    count = 0
    while frames:
        frame = frames[-1]
        if position == len(text):
            raise ValueError("Sequential bracket is not closed.")
        character = text[position]
        if character in "]}":
            closing = "]" if frame["opening"] == "[" else "}"
            if character != closing:
                raise ValueError(
                    "Sequential bracket has mismatched delimiter "
                    f"at offset {position}."
                )
            if not frame["entries"]:
                raise ValueError("Sequential bracket contains no valid entries.")
            entries = tuple(frame["entries"])
            if frame["legacy"]:
                entries = _legacy_entries(entries)
            scope = BracketArm(
                text[frame["start"] : position + 1],
                entries,
                frame["opening"],
                closing,
                normalized_legacy=frame["legacy"],
            )
            frames.pop()
            position += 1
            if not frames:
                if position != len(text):
                    raise ValueError(
                        "Sequential bracket has unrecognised trailing content."
                    )
                return scope
            frames[-1]["entries"].append(scope)
            continue

        start = position
        dotted = character == "."
        if dotted:
            position += 1
            if position == len(text):
                raise ValueError("Sequential bracket ends with an empty entry.")
            character = text[position]
        if character in "[{":
            if not dotted:
                if not frame["entries"]:
                    frame["legacy"] = True
                elif text[position + 1 : position + 2] != ".":
                    raise ValueError(
                        "Sequential bracket child requires .[...] or [.arm]."
                    )
            if len(frames) >= MAX_BRACKET_DEPTH:
                raise ValueError(
                    f"Sequential bracket nesting exceeds {MAX_BRACKET_DEPTH}; "
                    "use percent segments to reduce nesting."
                )
            frames.append(
                {
                    "start": start,
                    "opening": str(character),
                    "entries": [],
                    "legacy": False,
                }
            )
            position += 1
            continue
        if frame["entries"] and not dotted:
            raise ValueError(
                f"Sequential bracket has unrecognised content at offset {position}; "
                "expected .Entry(r,r) or a child scope."
            )
        match = _SCOPED_ENTRY_RE.match(text, position)
        if match is None or (match.group(2) is None and not match[1].startswith("!")):
            reason = (
                "has unrecognised content"
                if frame["entries"]
                else "contains no valid entries; unrecognised content"
            )
            raise ValueError(
                f"Sequential bracket {reason} at offset {position}: "
                f"{str(text[position : position + 40])!r}; expected Entry(r,r)."
            )
        slots = (match.group(2), match.group(3))
        if slots[0] is not None:
            slots = tuple(str(int(slot)) for slot in slots)
            if "0" in slots:
                raise ValueError(
                    "Sequential bracket attachment slots must be positive."
                )
        frame["entries"].append(
            BracketEntry(
                text[start : match.end()],
                _source_group(match, 1),
                *slots,
            )
        )
        position = match.end()
        count += 1
        if count > MAX_BRACKET_ENTRIES:
            raise ValueError(
                f"Sequential bracket exceeds {MAX_BRACKET_ENTRIES:,} entries; "
                "split it into smaller documents."
            )


def _bracket_end(text, start):
    """Return the end of one balanced bracket, rejecting mixed delimiters."""
    text = str(text)
    stack = []
    index = start
    while index < len(text):
        character = text[index]
        if character == "<":
            end = text.find(">", index + 1)
            if end < 0:
                raise ValueError("Synthetic monomer token is not closed.")
            index = end + 1
            continue
        if character in "[{":
            stack.append("]" if character == "[" else "}")
            if len(stack) > MAX_BRACKET_DEPTH:
                raise ValueError(
                    f"Sequential bracket nesting exceeds {MAX_BRACKET_DEPTH}; "
                    "use percent segments to reduce nesting."
                )
        elif character in "]}":
            if not stack or character != stack.pop():
                raise ValueError(
                    f"Sequential bracket has mismatched delimiter at offset {index}."
                )
            if not stack:
                return index + 1
        index += 1
    raise ValueError(f"Sequential bracket at offset {start} is not closed.")


def bracket_regions(text):
    """Yield complete .[...] / .{...} spans, including any nested arms."""
    text = str(text)
    position = 0
    while position < len(text) - 1:
        if text[position] == "<":
            end = text.find(">", position + 1)
            if end < 0:
                raise ValueError("Synthetic monomer token is not closed.")
            position = end + 1
        elif text[position] == "." and text[position + 1] in "[{":
            end = _bracket_end(text, position + 1)
            yield position, end
            position = end
        else:
            position += 1


def crosslink_declarations(entries):
    """Yield declared slot pairs in bracket traversal order, including arms."""
    pending = list(reversed(entries))
    while pending:
        entry = pending.pop()
        if isinstance(entry, BracketArm):
            pending.extend(reversed(entry.entries))
        elif entry.token.startswith("!") and entry.previous_slot is not None:
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
        previous, own = str(int(previous)), str(int(own))
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
                    f"(expected inverse R{first_own}→partner R{first_previous})."
                )
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


def _normalized_scope(scope):
    """Render only legacy-containing subtrees, retaining token source slices."""
    rendered = {}
    pending = [(scope, False)]
    while pending:
        node, ready = pending.pop()
        children = [entry for entry in node.entries if isinstance(entry, BracketArm)]
        if not ready:
            pending.append((node, True))
            pending.extend((child, False) for child in reversed(children))
            continue
        changed = (
            node.legacy
            or node.normalized_legacy
            or any(rendered[id(child)][1] for child in children)
        )
        if not changed:
            rendered[id(node)] = (node.text, False)
            continue
        raw = str(node.text)
        if node.legacy:
            prefix = "[."
        else:
            prefix = ("." if raw.startswith(".") else "") + node.opening
            opening = 1 if raw.startswith(".") else 0
            if raw[opening + 1 : opening + 2] == ".":
                prefix += "."
        pieces = [prefix]
        for position, entry in enumerate(node.entries):
            if isinstance(entry, BracketArm):
                pieces.append(rendered[id(entry)][0])
            else:
                text = entry.text
                if position == 0 and text.startswith("."):
                    text = text[1:]
                elif position and not text.startswith("."):
                    text = "." + text
                pieces.append(text)
        pieces.append(node.closing)
        rendered[id(node)] = (_source_join("", pieces), True)
    return rendered[id(scope)][0]


def normalize_legacy_brackets(seg):
    """Normalize legacy hub scopes, including those inside modern children.

    Parsing uses the same recursive records as Sequence. This editing helper
    writes their explicit arms only when a legacy wrapper needs normalization;
    it never determines graph ownership by repeatedly rewriting substrings.
    """
    out, previous_end = [], 0
    for start, end in bracket_regions(seg):
        scope = parse_bracket_group(seg[start:end])
        normalized = _normalized_scope(scope)
        if normalized != scope.text:
            out.extend((seg[previous_end:start], normalized))
            previous_end = end
    if not out:
        return seg
    out.append(seg[previous_end:])
    return _source_join("", out)


def _flatten_one_nested(s):
    """Compatibility entry point for one complete legacy hub scope."""
    if not s.startswith("[[") or _bracket_end(s, 0) != len(s):
        raise ValueError(
            "Legacy nested bracket must be one complete [[...].Entry(r,r)] group."
        )
    return _normalized_scope(parse_bracket_group(s))


def split_outside(string, by_element, outside, keep_marker=True):
    """
    Splits a string by delimiter only if outside of a given delimiter

    :param string: string to be split
    :param by_element: delimiter(s) by which to be split
    :param outside: only split if outside of this
    :param keep_marker: if True keep the chunk marker, remove otherwise

    :return splitChains: split string as list
    """
    delimiters = set(by_element)
    if len(outside) == 1:
        outside += outside
    pairs = {outside[0]: outside[1]}
    if outside == '[]':
        pairs['{'] = '}'
    stack, pieces, result = [], [], []
    for position in range(len(string)):
        character = string[position]
        if stack and character == stack[-1]:
            stack.pop()
            if keep_marker:
                pieces.append(character)
        elif character in pairs:
            stack.append(pairs[character])
            if keep_marker:
                pieces.append(character)
        elif not stack and character in delimiters:
            result.append(_source_join('', pieces))
            pieces = []
        else:
            pieces.append(character)
    result.append(_source_join('', pieces))
    return result
