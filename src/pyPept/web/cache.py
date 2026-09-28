"""Bounded cache for rendered molecules."""

from __future__ import annotations

from collections import OrderedDict as _OD
from threading import RLock

from pyPept import monomer_store

_render_cache = _OD()
_RENDER_CACHE_MAX = 200
_cache_lock = RLock()


def library_version():
    path = monomer_store.library_path()
    stamp = path.stat()
    return (str(path), stamp.st_mtime_ns, stamp.st_size)


def _rc_get(key):
    with _cache_lock:
        if key in _render_cache:
            _render_cache.move_to_end(key)
            return _render_cache[key]
    return None


def _rc_put(key, value):
    with _cache_lock:
        _render_cache[key] = value
        _render_cache.move_to_end(key)
        if len(_render_cache) > _RENDER_CACHE_MAX:
            _render_cache.popitem(last=False)
