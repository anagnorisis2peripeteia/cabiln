"""Source locations carried through the existing CABILN lowering.

This module has no notation grammar. ``Sequence`` records occurrences while its
existing parser lowers synthetic tokens, branches and inline attachments. Plain
strings retain the usual fast path; source tracking is opt-in for editing.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True, order=True)
class Span:
    start: int
    end: int


class SourceError(ValueError):
    """An input failure with a repair hint and an optional original source span."""

    def __init__(self, message, *, hint, span=None):
        super().__init__(message)
        self.hint = hint
        self.span = span


@dataclass(frozen=True)
class Occurrence:
    """Original token and attachment ownership for one resolved occurrence.

    ``scope`` is its actual containing scope, or a legacy entry's real span.
    Flat continuations share that scope and name their preceding token in
    ``host``. ``terminal`` is false if inserting a flat entry here would redirect
    any following marker, continuation, or child. ``bracket``/``arm`` retain the
    outermost/nearest-child locations expected by existing editing clients.
    """

    token: Span
    entry: Span
    kind: Literal["explicit", "inline", "bracket"]
    segment: int
    bracket: Span | None = None
    arm: Span | None = None
    terminal: bool = True
    protected: bool = False
    bracketable: bool = False
    legacy_arm: bool = False
    host: Span | None = None
    scope: Span | None = None
    parent_scope: Span | None = None


@dataclass(frozen=True)
class SourceGroup:
    """A parsed scope and the monomer that hosts it, before ID projection."""

    span: Span
    host: Span
    parent: Span | None
    kind: str
    opening: str
    closing: str
    protected: bool
    legacy: bool = False


@dataclass(frozen=True)
class BondMarker:
    span: Span
    slot: int
    owner: Span
    kind: Literal["inline", "terminal", "bracket"]
    bracket: Span | None = None
    arm: Span | None = None
    scope: Span | None = None


class Tracker:
    def __init__(self, source):
        self.source = source
        self.entries = {}
        self.roots = []
        self.segment = 0
        self.labels = set()
        self.markers = []
        self.groups = []

    def root(self, text, segment):
        span = origin_span(text)
        if span is not None:
            self.roots.append((span, segment))

    def entry(
        self,
        token,
        entry,
        kind,
        bracket=None,
        arm=None,
        terminal=True,
        *,
        host=None,
        scope=None,
        parent_scope=None,
        protected=None,
        legacy_arm=None,
    ):
        from pyPept.notation import supports_bracket_token

        span = origin_span(token)
        if span is not None:
            arm_span = origin_span(arm)
            self.entries[span] = Occurrence(
                span,
                origin_span(entry),
                kind,
                self.segment,
                origin_span(bracket),
                arm_span,
                terminal,
                (
                    bool(bracket is not None and str(bracket).startswith(".{"))
                    if protected is None
                    else protected
                ),
                supports_bracket_token(self.source[span.start : span.end]),
                (
                    arm_span is not None
                    and not self.source[arm_span.start : arm_span.end].startswith("[")
                    if legacy_arm is None
                    else legacy_arm
                ),
                origin_span(host),
                origin_span(scope),
                origin_span(parent_scope),
            )

    def occurrence(self, token):
        from pyPept.notation import supports_bracket_token

        span = origin_span(token)
        if span is None:
            raise ValueError("Lowered monomer lost its source")
        if span in self.entries:
            return self.entries[span]
        region = next(
            (
                (r, s)
                for r, s in self.roots
                if r.start <= span.start and span.end <= r.end
            ),
            None,
        )
        if region is None:
            raise ValueError("Lowered monomer has no original source region")
        root, segment = region
        return Occurrence(
            span,
            root,
            "explicit",
            segment,
            bracketable=supports_bracket_token(self.source[span.start : span.end]),
        )


class SourceText(str):
    def __new__(cls, value, origins, tracker):
        obj = str.__new__(cls, value)
        obj.origins = tuple(origins)
        obj.tracker = tracker
        assert len(obj) == len(obj.origins)
        return obj

    @classmethod
    def original(cls, value):
        return cls(value, (Span(i, i + 1) for i in range(len(value))), Tracker(value))

    def __getitem__(self, key):
        if isinstance(key, int):
            key = key if key >= 0 else len(self) + key
            if key < 0 or key >= len(self):
                raise IndexError("string index out of range")
            key = slice(key, key + 1)
        return SourceText(str.__getitem__(self, key), self.origins[key], self.tracker)

    def __add__(self, other):
        return join("", [self, other])

    def __radd__(self, other):
        return join("", [other, self])

    def split(self, sep=None, maxsplit=-1):
        if not sep:
            raise ValueError("SourceText.split requires a non-empty separator")
        result, start = [], 0
        while maxsplit != 0:
            pos = str.find(self, sep, start)
            if pos < 0:
                break
            result.append(self[start:pos])
            start = pos + len(sep)
            maxsplit -= 1
        result.append(self[start:])
        return result

    def strip(self, chars=None):
        left = len(self) - len(str.lstrip(self, chars))
        end = len(str.rstrip(self, chars))
        return self[left : max(left, end)]


def join(sep, pieces):
    pieces = list(pieces)
    mapped = next((p for p in pieces if isinstance(p, SourceText)), None)
    if mapped is None:
        return sep.join(pieces)
    origins = []
    for i, p in enumerate(pieces):
        if i:
            origins.extend([None] * len(sep))
        if isinstance(p, SourceText) and p.tracker is not mapped.tracker:
            raise ValueError("Cannot combine text from different source revisions")
        origins.extend(p.origins if isinstance(p, SourceText) else [None] * len(p))
    return SourceText(sep.join(map(str, pieces)), origins, mapped.tracker)


def group(match, num=0):
    return match.string[match.start(num) : match.end(num)]


def sub(pattern, replacement, value):
    if not isinstance(value, SourceText):
        return re.sub(pattern, replacement, value)
    pieces, end = [], 0
    for match in re.finditer(pattern, value):
        pieces.append(value[end : match.start()])
        pieces.append(
            replacement(match) if callable(replacement) else match.expand(replacement)
        )
        end = match.end()
    pieces.append(value[end:])
    return join("", pieces)


def origin_span(text):
    if not isinstance(text, SourceText):
        return None
    if hasattr(text, "_origin_span"):
        return text._origin_span
    origins = [p for p in text.origins if p is not None]
    text._origin_span = (
        Span(min(p.start for p in origins), max(p.end for p in origins))
        if origins
        else None
    )
    return text._origin_span


def synthetic(match, symbol):
    value = group(match)
    if not isinstance(value, SourceText):
        return symbol
    return SourceText(symbol, [origin_span(value)] * len(symbol), value.tracker)


def record(token, entry, kind, bracket=None, arm=None, terminal=True, **ownership):
    if isinstance(token, SourceText):
        token.tracker.entry(token, entry, kind, bracket, arm, terminal, **ownership)


def scope_group(
    value,
    *,
    host,
    parent=None,
    kind="bracket",
    opening="[",
    closing="]",
    protected=False,
    legacy=False,
):
    if not isinstance(value, SourceText):
        return
    span, owner = origin_span(value), origin_span(host)
    if span is None or owner is None:
        raise ValueError("Sequential bracket lost its source group or host")
    value.tracker.groups.append(
        SourceGroup(
            span,
            owner,
            origin_span(parent),
            kind,
            opening,
            closing,
            protected,
            legacy,
        )
    )


def bond_marker(
    value, slot, *, owner=None, bracket=None, arm=None, terminal=False, scope=None
):
    """Record an endpoint where the existing lowering resolves its owner/slot."""
    if not isinstance(value, SourceText):
        return
    span = origin_span(value)
    owner_span = origin_span(owner)
    if owner_span is None:
        owner_span = next(
            (
                root
                for root, _ in value.tracker.roots
                if root.start <= span.start and span.end <= root.end
            ),
            None,
        )
    if owner_span is None:
        raise ValueError("Crosslink endpoint has no original source owner")
    kind = "terminal" if terminal else "bracket" if bracket is not None else "inline"
    value.tracker.markers.append(
        BondMarker(
            span,
            slot,
            owner_span,
            kind,
            origin_span(bracket),
            origin_span(arm),
            origin_span(scope),
        )
    )
