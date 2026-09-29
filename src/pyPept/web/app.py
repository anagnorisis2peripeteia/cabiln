"""Application setup and command-line entry point for CABILN."""

from __future__ import annotations

import argparse
import asyncio
from contextlib import asynccontextmanager
import ipaddress
import json
import os
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from . import builder, conversion, monomers, projects, rendering
from .execution import (
    ChemistryExecutor, ExecutionConfig, ExecutionMiddleware,
    RequestLoggingMiddleware, _internal_error, configure_logging,
    deployment_version,
)
from .readiness import check_readiness
from .security import configure_registration, require_registration

STATIC_DIR = Path(__file__).resolve().parent / "static"
SERVER_ID = deployment_version()


def create_app(*, allow_registration=None, execution_mode=None, observability=True):
    config = ExecutionConfig.from_env(execution_mode)
    if (
        os.environ.get("CABILN_ENV") == "production"
        and config.mode != "process" and execution_mode is None
    ):
        raise ValueError("Production requires CABILN_EXECUTION=process")
    if os.environ.get("CABILN_ENV") == "production":
        from rdkit import rdBase

        # Authorized parent-process registration also parses user structures.
        # Expected errors remain in HTTP responses, never native stderr logs.
        rdBase.DisableLog("rdApp.*")
    version = deployment_version()

    @asynccontextmanager
    async def lifespan(app):
        try:
            app.state.readiness = await asyncio.to_thread(check_readiness, STATIC_DIR)
            if config.mode == "process":
                app.state.executor = ChemistryExecutor(config)
                await app.state.executor.start()
        except Exception as exc:
            app.state.readiness = {"ready": False}
            _internal_error(exc)
        try:
            yield
        finally:
            if app.state.executor is not None:
                await app.state.executor.close()

    app = FastAPI(title="CABILN peptide builder", lifespan=lifespan)
    app.state.readiness = {"ready": False}
    app.state.executor = None
    app.state.execution_config = config
    app.state.release = version
    app.state.allow_registration = (
        os.environ.get("CABILN_ENABLE_REGISTRATION") == "1"
        if allow_registration is None
        else allow_registration
    )
    configure_registration(app, explicit_trusted_local=(allow_registration is True))
    if config.mode == "process":
        app.add_middleware(ExecutionMiddleware, state=app.state, config=config)
    app.add_middleware(GZipMiddleware, minimum_size=1000, compresslevel=4)
    if observability:
        configure_logging()
        app.add_middleware(RequestLoggingMiddleware, version=version)

    @app.exception_handler(RequestValidationError)
    async def invalid_request(request, exc):
        errors = exc.errors()
        first = errors[0]
        field = ".".join(str(part) for part in first["loc"][1:])
        return JSONResponse({"error": f"{field}: {first['msg']}"}, status_code=422)

    app.mount(
        "/static", StaticFiles(directory=STATIC_DIR, check_dir=False), name="static"
    )
    for router in (
        builder.router,
        conversion.router,
        monomers.router,
        projects.router,
        rendering.router,
    ):
        app.include_router(router)

    @app.get("/", include_in_schema=False)
    def index():
        return FileResponse(STATIC_DIR / "index.html")

    @app.get("/register", include_in_schema=False)
    def register_page(request: Request):
        require_registration(request)
        return FileResponse(STATIC_DIR / "register.html")

    @app.get("/capabilities")
    def capabilities():
        return {"registration": app.state.allow_registration}

    @app.get("/server_id", response_class=PlainTextResponse)
    def server_id():
        return version

    @app.get("/examples")
    def examples():
        return json.loads((STATIC_DIR / "examples.json").read_text(encoding="utf-8"))

    @app.get("/health")
    async def health():
        return {"status": "ok"}

    @app.get("/ready")
    async def ready():
        available = app.state.readiness.get("ready", False)
        if config.mode == "process":
            executor = app.state.executor
            available = available and executor is not None and executor.ready
        return JSONResponse(
            {"status": "ready" if available else "not_ready", "release": version},
            status_code=200 if available else 503,
        )

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
        configure_registration(app, explicit_trusted_local=True)
    try:
        loopback = ipaddress.ip_address(args.host).is_loopback
    except ValueError:
        loopback = args.host == "localhost"
    if (
        app.state.allow_registration and not app.state.registration_token
        and not loopback
    ):
        parser.error("Public registration requires CABILN_REGISTRATION_TOKEN")
    options = {}
    if app.state.execution_config.mode == "process":
        options["limit_concurrency"] = 16
        options["workers"] = 1
    uvicorn.run(app, host=args.host, port=args.port, access_log=False, **options)
