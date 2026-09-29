"""Bounded cache for rendered molecules."""

from __future__ import annotations

from collections import OrderedDict as _OD
import os
import sys
from threading import RLock

from pyPept import monomer_store

_render_cache = _OD()
_RENDER_CACHE_MAX = 200
_RENDER_CACHE_BYTES = int(os.environ.get("CABILN_CACHE_BYTES", 32 * 1024 * 1024))
if not 0 <= _RENDER_CACHE_BYTES <= 512 * 1024 * 1024:
    raise ValueError("CABILN_CACHE_BYTES must be between 0 and 512 MiB")
_cache_lock = RLock()
_entry_sizes = {}
_cache_bytes = 0


def _retained_size(value):
    """Account for Python payload containers, strings and scalars, once each."""
    pending, seen, total = [value], set(), 0
    while pending:
        item = pending.pop()
        if id(item) in seen:
            continue
        seen.add(id(item))
        total += sys.getsizeof(item)
        if isinstance(item, dict):
            pending.extend(item.keys())
            pending.extend(item.values())
        elif isinstance(item, (list, tuple)):
            pending.extend(item)
    return total


def clear_render_cache():
    global _cache_bytes
    with _cache_lock:
        _render_cache.clear()
        _entry_sizes.clear()
        _cache_bytes = 0


def library_version():
    return monomer_store.library_version()


def _rc_get(key):
    with _cache_lock:
        if key in _render_cache:
            _render_cache.move_to_end(key)
            return _render_cache[key]
    return None


def _rc_put(key, value):
    global _cache_bytes
    size = _retained_size((key, value))
    with _cache_lock:
        if key in _render_cache:
            del _render_cache[key]
            _cache_bytes -= _entry_sizes.pop(key)
        if size > _RENDER_CACHE_BYTES:
            return
        _render_cache[key] = value
        _entry_sizes[key] = size
        _cache_bytes += size
        _render_cache.move_to_end(key)
        while (
            len(_render_cache) > _RENDER_CACHE_MAX
            or _cache_bytes > _RENDER_CACHE_BYTES
        ):
            oldest, _ = _render_cache.popitem(last=False)
            _cache_bytes -= _entry_sizes.pop(oldest)
