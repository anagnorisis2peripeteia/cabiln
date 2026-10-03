"""Private, bounded library snapshots addressed by opaque tab tokens."""

from __future__ import annotations

import secrets
import shutil
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from fastapi import HTTPException, Request

from pyPept.monomer_store import (
    _load_sdf, library_path, library_version, monomer_table, use_library,
)

from .execution import CHEMISTRY_ROUTES, error_response
from .security import require_same_origin

SESSION_HEADER = "X-Cabiln-Library"


@dataclass
class _Snapshot:
    path: Path
    base: tuple
    touched: float
    active: int = 1
    ready: bool = False


class SessionLibraries:
    """Keep at most 32 snapshots, expiring after an hour without a request.

    Active requests pin their snapshot. Browser-held definitions can recreate an
    evicted snapshot. Worker renewal retains the parent-owned files.
    """

    def __init__(self, maximum=32, lifetime=3600):
        self.maximum, self.lifetime = maximum, lifetime
        self.entries = {}
        self.directory = None

    def close(self):
        if self.directory is not None:
            self.directory.cleanup()
            self.directory = None
        self.entries.clear()

    def _remove(self, token):
        entry = self.entries.pop(token)
        shutil.rmtree(entry.path.parent)

    def _expire(self):
        now = time.monotonic()
        for token, entry in list(self.entries.items()):
            if not entry.active and now - entry.touched > self.lifetime:
                self._remove(token)

    def create(self):
        self._expire()
        if len(self.entries) >= self.maximum:
            idle = [(entry.touched, token) for token, entry in self.entries.items()
                    if not entry.active]
            if not idle:
                raise HTTPException(503, "Temporary libraries are busy; retry shortly.",
                                    headers={"Retry-After": "1"})
            self._remove(min(idle)[1])
        if self.directory is None:
            self.directory = tempfile.TemporaryDirectory(prefix="cabiln-libraries-")
        token = secrets.token_hex(24)
        directory = Path(self.directory.name) / token
        directory.mkdir(mode=0o700)
        entry = _Snapshot(directory / "monomers.sdf", library_version(),
                          time.monotonic())
        self.entries[token] = entry
        return token, entry

    def publish(self, token):
        entry = self.entries[token]
        if not entry.path.is_file() or entry.base != library_version():
            raise HTTPException(409, "The shared library changed; add the monomer again.")
        entry.ready = True

    def release(self, token):
        entry = self.entries[token]
        entry.active -= 1
        entry.touched = time.monotonic()
        if not entry.ready and not entry.active:
            self._remove(token)

    @contextmanager
    def select(self, token):
        self._expire()
        entry = self.entries.get(token)
        if entry is None or not entry.ready or entry.base != library_version():
            raise HTTPException(
                409, "This temporary library expired; restore its tab definitions.",
                headers={"X-Cabiln-Library-Expired": "1"},
            )
        entry.active += 1
        try:
            with use_library(entry.path):
                yield entry.path
        finally:
            self.release(token)


def write_overlay(destination, molecules):
    """Validate collisions once and write a private copy of the SDF/alias pair."""
    from rdkit import Chem
    from pyPept.library_quality import definition_hash

    source, version = library_path(), library_version()
    existing, by_name = _load_sdf()
    table = monomer_table()
    reserved = set(by_name)
    for key in ("_synonyms", "_degen_aliases", "_ambiguous_aliases"):
        reserved.update(table.attrs.get(key, {}))
    additions = []
    if sum(mol.GetNumAtoms() for mol in molecules) > 4096:
        raise ValueError("Custom monomers in one tab cannot exceed 4,096 atoms in total")
    names = set()
    for mol in molecules:
        name = mol.GetProp("symbol")
        if name in names:
            raise ValueError(f"Duplicate temporary monomer '{name}'")
        names.add(name)
        if name in reserved:
            if name in by_name and definition_hash(mol) == definition_hash(by_name[name]):
                continue
            raise ValueError(f"Monomer name '{name}' already exists with another definition or as an alias")
        additions.append(mol)
    destination = Path(destination)
    shutil.copyfile(source, destination)
    aliases = source.with_name("monomers.csv")
    if aliases.is_file():
        shutil.copyfile(aliases, destination.with_name("monomers.csv"))
    with destination.open("ab") as handle:
        with source.open("rb") as original:
            if source.stat().st_size:
                original.seek(-1, 2)
                if original.read(1) != b"\n":
                    handle.write(b"\n")
        for mol in additions:
            handle.write(Chem.SDWriter.GetText(mol, kekulize=False).encode("utf-8"))
    if library_version() != version:
        raise ValueError("The shared library changed; add the monomer again")
    return len(existing) + len(additions)


class SessionLibraryMiddleware:
    def __init__(self, app, libraries):
        self.app, self.libraries = app, libraries

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        route = scope["method"], scope["path"]
        if route not in CHEMISTRY_ROUTES and route != ("GET", "/project_context"):
            await self.app(scope, receive, send)
            return
        request = Request(scope)

        async def private_send(event):
            if event["type"] == "http.response.start":
                event["headers"] = [*event.get("headers", []),
                                    (b"cache-control", b"no-store")]
            await send(event)

        try:
            if route == ("POST", "/session_library"):
                require_same_origin(request)
                token, entry = self.libraries.create()
                scope["cabiln.session_destination"] = (token, str(entry.path))

                async def publish_send(event):
                    if event["type"] == "http.response.start" and event["status"] < 300:
                        self.libraries.publish(token)
                    await private_send(event)

                try:
                    await self.app(scope, receive, publish_send)
                finally:
                    self.libraries.release(token)
            elif token := request.headers.get(SESSION_HEADER):
                with self.libraries.select(token) as path:
                    scope["cabiln.library_path"] = str(path)
                    await self.app(scope, receive, private_send)
            else:
                await self.app(scope, receive, send)
        except HTTPException as exc:
            await error_response(exc)(scope, receive, private_send)
