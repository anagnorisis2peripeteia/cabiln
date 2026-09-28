"""Notation normalization shared by builder and conversion requests."""

from __future__ import annotations

import re


def split_top_level(text: str, delimiter: str) -> list[str]:
    """Split notation while preserving brackets and inline SMILES tokens."""
    parts, stack, start = [], [], 0
    closing = {"[": "]", "{": "}", "(": ")", "<": ">"}
    for index, char in enumerate(text):
        if stack and stack[-1] == ">":
            if char == ">":
                stack.pop()
        elif char in closing:
            stack.append(closing[char])
        elif char in closing.values():
            if not stack or stack.pop() != char:
                raise ValueError("Unbalanced brackets in CABILN")
        elif char == delimiter and not stack:
            parts.append(text[start:index])
            start = index + 1
    if stack:
        raise ValueError("Unbalanced brackets in CABILN")
    parts.append(text[start:])
    return parts


def backbone_token_indices(tokens: list[str]) -> list[int]:
    """Terminal !n markers are bonds, and do not consume a residue index."""
    return [
        i
        for i, token in enumerate(tokens)
        if token and not re.fullmatch(r"![A-Za-z0-9_]+", token)
    ]


def _to_bracket(cabiln: str) -> str:
    """Normalise positional-% branch notation to bracket form for Sequence.

    Pure crosslink branches (%TBMB.!1.!2.!3, %C20FA-AEEA-E_g.!1, etc.) are
    already parseable by Sequence and are left untouched.  Only branches where
    every segment lacks an !n marker need conversion.
    """
    if "%" not in cabiln:
        return cabiln
    branch_parts = split_top_level(cabiln, "%")[1:]
    if all("!" in p for p in branch_parts):
        return cabiln  # all crosslink — Sequence handles it
    from pyPept.sequence import cabiln_to_bracket

    return cabiln_to_bracket(cabiln)


def _renumber_xlinks(cabiln: str) -> str:
    """Renumber !X tags so they appear as !1, !2, … in left-to-right reading order."""
    seen: dict[int, int] = {}

    def replace(match):
        old = int(match.group(1))
        if old not in seen:
            seen[old] = len(seen) + 1
        return f"!{seen[old]}"

    return re.sub(r"!(\d+)(?![A-Za-z0-9_])", replace, cabiln)


def _apply_notation(cabiln: str, notation: str) -> str:
    if notation == "bracket":
        from pyPept.sequence import cabiln_to_bracket

        return _verified_notation(cabiln, _renumber_xlinks(cabiln_to_bracket(cabiln)))
    return cabiln


def parse_source(cabiln: str, warning_sink=None):
    """Return the parsed source explicitly, including legacy positional input.

    Callers exposing residue selections must use the returned source for edits;
    legacy conversion can change instance ordering.
    """
    from pyPept.sequence import Sequence

    messages = []
    parsed_source = cabiln
    try:
        sequence = Sequence(cabiln, warning_sink=messages.append)
    except ValueError:
        parsed_source = _to_bracket(cabiln)
        if parsed_source == cabiln:
            raise
        messages.clear()
        sequence = Sequence(parsed_source, warning_sink=messages.append)
        messages.insert(0, "Converted legacy positional notation to bracket form.")
    if warning_sink is not None:
        for message in messages:
            warning_sink(message)
    return sequence, parsed_source


def _verified_notation(source: str, result: str) -> str:
    """A notation conversion must retain the exact assembled structure."""
    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence
    from pyPept.structure import compare_structures

    sequence, _ = parse_source(source)
    original = Molecule(sequence).get_molecule(fmt="ROMol")
    if source != result:
        converted = Molecule(
            Sequence(result, warning_sink=lambda message: None)
        ).get_molecule(fmt="ROMol")
        if not compare_structures(original, converted).exact:
            raise ValueError("Notation conversion would change the molecular structure")
    return result
