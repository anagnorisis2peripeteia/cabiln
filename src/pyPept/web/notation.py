"""Compatibility imports for notation workflows now owned by the core package."""

from __future__ import annotations

import re

from pyPept.inputs import _to_bracket, format_source, parse_source, split_top_level


def backbone_token_indices(tokens: list[str]) -> list[int]:
    """Terminal !n markers are bonds, and do not consume a residue index."""
    return [
        i
        for i, token in enumerate(tokens)
        if token and not re.fullmatch(r"![A-Za-z0-9_]+", token)
    ]
