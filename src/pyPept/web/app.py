"""Application setup and command-line entry point for CABILN."""

from __future__ import annotations

import argparse
import json
import os
import uuid
from pathlib import Path

from fastapi import FastAPI
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from . import builder, conversion, monomers, rendering

STATIC_DIR = Path(__file__).resolve().parent / "static"
SERVER_ID = uuid.uuid4().hex[:8]


def create_app(*, allow_registration=None):
    app = FastAPI(title="CABILN peptide builder")
    app.state.allow_registration = (
        os.environ.get("CABILN_ENABLE_REGISTRATION") == "1"
        if allow_registration is None
        else allow_registration
    )

    @app.exception_handler(RequestValidationError)
    async def invalid_request(request, exc):
        errors = exc.errors()
        first = errors[0]
        field = ".".join(str(part) for part in first["loc"][1:])
        return JSONResponse({"error": f"{field}: {first['msg']}"}, status_code=422)

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    for router in (
        builder.router,
        conversion.router,
        monomers.router,
        rendering.router,
    ):
        app.include_router(router)

    @app.get("/", include_in_schema=False)
    def index():
        return FileResponse(STATIC_DIR / "index.html")

    @app.get("/register", include_in_schema=False)
    def register_page():
        if not app.state.allow_registration:
            return PlainTextResponse(
                "This monomer library is read-only.", status_code=403
            )
        return FileResponse(STATIC_DIR / "register.html")

    @app.get("/capabilities")
    def capabilities():
        return {"registration": app.state.allow_registration}

    @app.get("/server_id", response_class=PlainTextResponse)
    def server_id():
        return SERVER_ID

    @app.get("/examples")
    def examples():
        return json.loads((STATIC_DIR / "examples.json").read_text(encoding="utf-8"))

    @app.get("/health")
    def health():
        return {"status": "ok"}

    return app


app = create_app()


def main():
    """Run locally; Render supplies PORT and an explicit public bind address."""
    import uvicorn

    parser = argparse.ArgumentParser(description="Run the CABILN peptide builder")
    parser.add_argument(
        "--host",
        default=os.environ.get(
            "HOST", "0.0.0.0" if "PORT" in os.environ else "127.0.0.1"
        ),
    )
    parser.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8732")))
    parser.add_argument(
        "--enable-registration",
        action="store_true",
        help="Allow library changes on a trusted local instance",
    )
    args = parser.parse_args()
    if args.enable_registration:
        app.state.allow_registration = True
    uvicorn.run(app, host=args.host, port=args.port)
