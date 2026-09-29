"""Access policy for the shared library's administrative writer."""

from __future__ import annotations

import base64
import binascii
import hmac
import ipaddress
import os
from pathlib import Path
from urllib.parse import urlsplit

from fastapi import HTTPException, Request


def configure_registration(app, *, explicit_trusted_local=False):
    """Keep public registration authenticated and local administration explicit."""
    token = os.environ.get("CABILN_REGISTRATION_TOKEN", "")
    production = os.environ.get("CABILN_ENV") == "production"
    app.state.registration_token = token
    app.state.registration_user = os.environ.get("CABILN_REGISTRATION_USER", "admin")
    app.state.registration_trusted_local = explicit_trusted_local
    if not app.state.allow_registration:
        return
    if token and not 32 <= len(token) <= 4096:
        raise ValueError("CABILN_REGISTRATION_TOKEN must contain 32–4096 characters")
    if production:
        if not token:
            raise ValueError("Production registration requires an administrative token")
        library = os.environ.get("CABILN_MONOMER_LIBRARY")
        backup = os.environ.get("CABILN_LIBRARY_BACKUP_DIR")
        if not library or not backup:
            raise ValueError(
                "Production registration requires an external library and backup "
                "directory (CABILN_MONOMER_LIBRARY, CABILN_LIBRARY_BACKUP_DIR)"
            )
        package = Path(__file__).resolve().parents[1]
        for path in (Path(library).resolve(), Path(backup).resolve()):
            if path == package or package in path.parents:
                raise ValueError("Administrative data must be outside the installation")


def _is_loopback(request):
    host = request.client.host if request.client else ""
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return host == "testclient" and request.app.state.registration_trusted_local


def require_registration(request: Request):
    """Authenticate each write; a page visit does not grant write permission."""
    state = request.app.state
    if not state.allow_registration:
        raise HTTPException(403, "This monomer library is read-only.")
    token = state.registration_token
    if token:
        credentials = request.headers.get("authorization", "")
        user, password = "", ""
        try:
            scheme, encoded = credentials.split(" ", 1)
            if scheme.lower() == "basic":
                decoded = base64.b64decode(encoded, validate=True).decode("utf-8")
                user, password = decoded.split(":", 1)
        except (ValueError, UnicodeError, binascii.Error):
            pass
        user_ok = hmac.compare_digest(
            user.encode("utf-8"), state.registration_user.encode("utf-8")
        )
        password_ok = hmac.compare_digest(password.encode("utf-8"), token.encode())
        if not (user_ok and password_ok):
            raise HTTPException(
                401,
                "Administrative authentication is required.",
                headers={"WWW-Authenticate": 'Basic realm="CABILN administration"'},
            )
    elif not _is_loopback(request):
        raise HTTPException(403, "Unauthenticated registration is local-only.")
    if request.method not in {"GET", "HEAD"}:
        origin = request.headers.get("origin")
        if origin:
            try:
                parsed = urlsplit(origin)
            except ValueError as exc:
                raise HTTPException(403, "Invalid registration origin.") from exc
            if (
                parsed.scheme not in {"http", "https"}
                or parsed.netloc.lower() != request.headers.get("host", "").lower()
            ):
                raise HTTPException(403, "Cross-origin registration is not permitted.")
        if request.headers.get("sec-fetch-site") == "cross-site":
            raise HTTPException(403, "Cross-origin registration is not permitted.")
