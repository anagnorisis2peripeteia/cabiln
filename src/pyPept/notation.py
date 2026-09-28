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
    i = 0
    while i < len(seg):
        if seg[i : i + 2] == ".[" and i + 2 < len(seg) and seg[i + 2] == "[":
            depth = 0
            j = i + 1
            while j < len(seg):
                if seg[j] == "[":
                    depth += 1
                elif seg[j] == "]":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            whole = seg[i + 1 : j + 1]
            flat = _flatten_one_nested(whole)
            out.append(seg[i : i + 1] + flat)
            i = j + 1
        else:
            out.append(seg[i])
            i += 1
    return _source_join("", out)


_ENTRY_ANY = re.compile("(" + _ENTRY + ")")


def _flatten_one_nested(s):
    """Convert a single nested bracket to new sub-bracket arm notation.

    Input ``s`` is the bracket expression starting at the outer ``[``.
    Returns ``[anchor[.inner_rest...][.outer...]]`` — every arm from the
    anchor is an explicit ``[.Entry(r,r)]`` sub-bracket so the parser's
    pointer semantics are unambiguous.
    """
    inner_start = s.index("[", 1)
    depth = 0
    inner_end = inner_start
    while inner_end < len(s):
        if s[inner_end] == "[":
            depth += 1
        elif s[inner_end] == "]":
            depth -= 1
            if depth == 0:
                break
        inner_end += 1

    inner_bracket = s[inner_start : inner_end + 1]
    inner_content = inner_bracket[1:-1]
    outer_part = s[inner_end + 1 : -1]
    outer_entries = [_source_group(m, 1) for m in _ENTRY_ANY.finditer(outer_part)]

    # If inner content already uses sub-bracket notation, preserve it and just
    # append outer entries as additional arms.
    if "[" in inner_content:
        outer_arms = _source_join("", ("[." + e + "]" for e in outer_entries))
        return s[:1] + inner_content + outer_arms + s[-1:]

    # Legacy flat inner content: split anchor from remaining chain entries.
    inner_entries = [_source_group(m, 1) for m in _ENTRY_ANY.finditer(inner_content)]
    anchor = inner_entries[0] if inner_entries else ""
    inner_rest = inner_entries[1:]

    inner_arms = _source_join("", ("[." + e + "]" for e in inner_rest))
    outer_arms = _source_join("", ("[." + e + "]" for e in outer_entries))
    return s[:1] + anchor + inner_arms + outer_arms + s[-1:]
