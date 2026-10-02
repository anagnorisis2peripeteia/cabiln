"""Bounded production chemistry jobs and request diagnostics.

Only an explicit set of read/edit HTTP operations crosses this process seam.
Workers execute the existing ASGI application, preserving its validation. There
is no waiting job queue: admission, IPC sizes, deadlines and worker replacement
belong to the parent. Static assets, health and administrator writes stay local.
"""

from __future__ import annotations

import asyncio
import base64
import contextvars
import json
import logging
import multiprocessing
import os
import socket
import struct
import sys
import time
import traceback
import uuid
from dataclasses import dataclass

from fastapi import HTTPException, Request
from fastapi.responses import JSONResponse, Response


CHEMISTRY_ROUTES = frozenset({
    ("POST", path) for path in (
        "/render", "/render_smiles", "/render_reference", "/render_mol",
        "/verify", "/convert_notation", "/smiles_to_cabiln", "/to_cabiln",
        "/preview_monomer", "/insert_bond", "/insert_backbone", "/validate_bond",
        "/replacement_options", "/replace_monomer",
        "/prepare_project", "/validate_project",
    )
} | {("GET", path) for path in (
    "/monomers", "/monomer_svg", "/monomer_rgroups", "/reactions",
)})
_request_id = contextvars.ContextVar("cabiln_request_id", default=None)
_diagnostics = contextvars.ContextVar("cabiln_diagnostics", default=None)
_logger = logging.getLogger("cabiln.runtime")


def deployment_version():
    from pyPept import __version__

    value = os.environ.get("CABILN_RELEASE") or os.environ.get("RENDER_GIT_COMMIT")
    return (value or __version__)[:128]


def _log(event, **fields):
    record = {"event": event, "request_id": _request_id.get(), **fields}
    _logger.info(json.dumps(record, sort_keys=True, separators=(",", ":")))


def configure_logging():
    if not _logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        _logger.addHandler(handler)
    _logger.setLevel(logging.INFO)


def _internal_error(exc):
    # Exception messages and source excerpts can contain submitted structures.
    # Keep only exception type and code locations, never locals, bodies or query.
    diagnostic = {
        "error_type": type(exc).__name__,
        "frames": [
            {"file": os.path.basename(frame.filename), "line": frame.lineno,
             "function": frame.name}
            for frame in traceback.extract_tb(exc.__traceback__)[-12:]
        ],
    }
    collected = _diagnostics.get()
    if collected is not None:
        if len(collected) < 8:
            collected.append(diagnostic)
    else:
        _log("internal_error", **diagnostic)


def error_response(exc):
    """Preserve expected input failures and sanitize unexpected server failures."""
    if isinstance(exc, HTTPException):
        return JSONResponse({"error": str(exc.detail)}, status_code=exc.status_code,
                            headers=exc.headers)
    if isinstance(exc, ValueError):
        from pyPept.source import SourceError

        message = str(exc).split("\n")[0]
        if message.startswith("Old BILN crosslink notation detected:"):
            message = (
                "This input uses legacy BILN crosslinks. Select BILN as the "
                "input notation, then use Convert to CABILN."
            )
        payload = {"error": message}
        if isinstance(exc, SourceError):
            payload["hint"] = exc.hint
            if exc.span is not None:
                payload["source_span"] = {"start": exc.span.start, "end": exc.span.end}
        return JSONResponse(payload, status_code=400)
    _internal_error(exc)
    return JSONResponse(
        {"error": "An internal error occurred.", "request_id": _request_id.get()},
        status_code=500,
    )


@dataclass(frozen=True)
class ExecutionConfig:
    mode: str = "local"
    workers: int = 1
    timeout_seconds: float = 30.0
    startup_seconds: float = 30.0
    memory_mb: int = 0
    max_request_bytes: int = 2 * 1024 * 1024
    max_response_bytes: int = 8 * 1024 * 1024
    max_query_bytes: int = 65536
    max_jobs: int = 100

    def __post_init__(self):
        if self.mode not in ("local", "process"):
            raise ValueError("CABILN_EXECUTION must be local or process")
        if not 1 <= self.workers <= 4:
            raise ValueError("CABILN_WORKERS must be between 1 and 4")
        if not 0 < self.timeout_seconds <= 120 or self.startup_seconds <= 0:
            raise ValueError("Job deadlines must be positive and at most 120 seconds")
        if not 0 <= self.memory_mb <= 4096:
            raise ValueError("Worker memory limit must be between 0 and 4096 MiB")
        if not 1024 <= self.max_request_bytes <= 4 * 1024 * 1024:
            raise ValueError("Request limit must be between 1 KiB and 4 MiB")
        if not 1024 <= self.max_response_bytes <= 32 * 1024 * 1024:
            raise ValueError("Response limit must be between 1 KiB and 32 MiB")
        if not 1024 <= self.max_query_bytes <= 65536:
            raise ValueError("Query limit must be between 1 and 64 KiB")
        if self.max_jobs < 1:
            raise ValueError("Worker recycling count must be positive")

    @classmethod
    def from_env(cls, mode=None):
        production = os.environ.get("CABILN_ENV") == "production"
        selected = mode or os.environ.get("CABILN_EXECUTION", "local")
        memory = int(os.environ.get(
            "CABILN_WORKER_MEMORY_MB", "1024" if production else "0"
        ))
        if production and selected == "process" and not memory:
            raise ValueError("Production requires a positive worker memory limit")
        return cls(
            mode=selected,
            workers=int(os.environ.get("CABILN_WORKERS", "1")),
            timeout_seconds=float(os.environ.get("CABILN_JOB_TIMEOUT_SECONDS", "30")),
            memory_mb=memory,
            max_request_bytes=int(os.environ.get(
                "CABILN_MAX_REQUEST_BYTES", str(2 * 1024 * 1024)
            )),
            max_response_bytes=int(os.environ.get(
                "CABILN_MAX_RESPONSE_BYTES", str(8 * 1024 * 1024)
            )),
        )


class ExecutionFailure(Exception):
    def __init__(self, message, status=503):
        super().__init__(message)
        self.status = status


def _encode(value):
    return json.dumps(value, separators=(",", ":")).encode("utf-8")


def _read_exact(sock, size):
    chunks = bytearray()
    while len(chunks) < size:
        part = sock.recv(size - len(chunks))
        if not part:
            raise EOFError("Worker transport closed")
        chunks.extend(part)
    return bytes(chunks)


def _read_message(sock, maximum):
    size = struct.unpack("!I", _read_exact(sock, 4))[0]
    if size > maximum:
        raise ValueError("Worker transport message exceeds its limit")
    return json.loads(_read_exact(sock, size))


def _send_message(sock, value):
    data = _encode(value)
    sock.sendall(struct.pack("!I", len(data)) + data)


async def _read_async(sock, maximum):
    loop = asyncio.get_running_loop()

    async def exact(size):
        chunks = bytearray()
        while len(chunks) < size:
            part = await loop.sock_recv(sock, size - len(chunks))
            if not part:
                raise EOFError("Worker transport closed")
            chunks.extend(part)
        return bytes(chunks)

    size = struct.unpack("!I", await exact(4))[0]
    if size > maximum:
        raise ValueError("Worker transport message exceeds its limit")
    return json.loads(await exact(size))


async def _serve_worker_request(app, message, config):
    body = base64.b64decode(message["body"], validate=True)
    query = base64.b64decode(message["query"], validate=True)
    if len(body) > config.max_request_bytes or len(query) > config.max_query_bytes:
        raise ValueError("Worker request exceeds its limit")
    if (message["method"], message["path"]) not in CHEMISTRY_ROUTES:
        raise ValueError("Operation is not permitted in a chemistry worker")
    scope = {
        "type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1",
        "method": message["method"], "scheme": "http", "path": message["path"],
        "raw_path": message["path"].encode(), "query_string": query,
        "root_path": "", "headers": [
            (b"content-type", message["content_type"].encode("latin1"))
        ],
        "client": ("127.0.0.1", 0), "server": ("worker", 0),
    }
    sent = False
    output = bytearray()
    status, headers = 500, []

    async def receive():
        nonlocal sent
        if not sent:
            sent = True
            return {"type": "http.request", "body": body, "more_body": False}
        await asyncio.Future()

    async def send(event):
        nonlocal status, headers
        if event["type"] == "http.response.start":
            status = event["status"]
            headers = event.get("headers", [])
        elif event["type"] == "http.response.body":
            part = event.get("body", b"")
            if len(output) + len(part) > config.max_response_bytes:
                raise ExecutionFailure("Chemistry response exceeds its byte limit", 413)
            output.extend(part)

    token = _request_id.set(message["request_id"])
    diagnostics = []
    diagnostic_token = _diagnostics.set(diagnostics)
    try:
        try:
            await app(scope, receive, send)
        except ExecutionFailure as exc:
            response = JSONResponse({"error": str(exc)}, status_code=exc.status)
            status, headers = response.status_code, response.raw_headers
            output = bytearray(response.body)
        except Exception as exc:
            response = error_response(exc)
            status, headers = response.status_code, response.raw_headers
            output = bytearray(response.body)
        return {
            "status": status,
            "headers": [(key.decode("latin1"), value.decode("latin1"))
                        for key, value in headers],
            "body": base64.b64encode(output).decode("ascii"),
            "diagnostics": diagnostics,
        }
    finally:
        _request_id.reset(token)
        _diagnostics.reset(diagnostic_token)


def _worker_main(sock, config):
    # Native chemistry warnings can echo user structures. Protocol data uses a
    # separate socket; incidental stdout/stderr never enters routine server logs.
    with open(os.devnull, "wb") as sink:
        os.dup2(sink.fileno(), 1)
        os.dup2(sink.fileno(), 2)
    try:
        from .app import STATIC_DIR, create_app
        from .monomers import palette_response
        from .readiness import check_readiness

        app = create_app(allow_registration=False, execution_mode="local",
                         observability=False)
        if config.memory_mb:
            if not sys.platform.startswith("linux"):
                raise RuntimeError("Worker memory enforcement requires Linux")
            import resource

            cap = config.memory_mb * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        check_readiness(STATIC_DIR)
        palette_response()
        _send_message(sock, {"ready": True})
        runner_class = getattr(asyncio, "Runner", _EventLoop)
        with runner_class() as runner:
            while True:
                message = _read_message(
                    sock,
                    config.max_request_bytes * 2 + config.max_query_bytes * 2 + 4096,
                )
                response = runner.run(_serve_worker_request(app, message, config))
                _send_message(sock, response)
    except (EOFError, BrokenPipeError, ConnectionResetError):
        pass
    except Exception as exc:
        try:
            _send_message(sock, {"ready": False, "error_type": type(exc).__name__})
        except OSError:
            pass
    finally:
        sock.close()


class _EventLoop:
    """The small Runner subset required on Python 3.9 and 3.10."""
    def __enter__(self):
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        return self

    def run(self, coroutine):
        return self.loop.run_until_complete(coroutine)

    def __exit__(self, *args):
        self.loop.run_until_complete(self.loop.shutdown_asyncgens())
        self.loop.run_until_complete(self.loop.shutdown_default_executor())
        self.loop.close()


@dataclass
class _Worker:
    process: object = None
    socket: object = None
    ready: bool = False
    busy: bool = False
    jobs: int = 0


class ChemistryExecutor:
    """An event-loop-owned bounded set of persistent, replaceable processes."""
    def __init__(self, config):
        self.config = config
        self.workers = [_Worker() for _ in range(config.workers)]
        self.closing = False
        self.tasks = set()
        self.recoveries = set()

    @property
    def ready(self):
        return not self.closing and any(
            worker.ready and worker.process is not None and worker.process.is_alive()
            for worker in self.workers
        )

    @property
    def available(self):
        return not self.closing and any(
            worker.ready and not worker.busy for worker in self.workers
        )

    async def start(self):
        # Set these before spawn imports NumPy/RDKit. The parent already imported
        # chemistry, but children must not create a BLAS thread pool per job slot.
        for variable in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
        ):
            os.environ[variable] = "1"
        starting = [asyncio.create_task(self._start(worker)) for worker in self.workers]
        try:
            await asyncio.gather(*starting)
        except BaseException:
            for task in starting:
                task.cancel()
            await asyncio.gather(*starting, return_exceptions=True)
            await self.close()
            raise

    async def _start(self, worker):
        parent, child = socket.socketpair()
        parent.setblocking(False)
        process = multiprocessing.get_context("spawn").Process(
            target=_worker_main, args=(child, self.config), daemon=True,
        )
        worker.process, worker.socket = process, parent
        try:
            process.start()
        finally:
            child.close()
        message = await asyncio.wait_for(
            _read_async(parent, 4096), timeout=self.config.startup_seconds
        )
        if message != {"ready": True}:
            _log("worker_startup_failed", error_type=message.get("error_type"))
            raise RuntimeError("Chemistry worker startup failed")
        worker.ready, worker.jobs = True, 0

    async def _stop(self, worker):
        worker.ready = False
        if worker.socket is not None:
            worker.socket.close()
            worker.socket = None
        process = worker.process
        if process is not None and process.pid is not None:
            if process.is_alive():
                process.kill()
            await asyncio.to_thread(process.join, 2)
            if process.is_alive():
                raise RuntimeError("Chemistry worker did not terminate")
            process.close()
        worker.process = None

    async def _retire(self, worker):
        # Termination completes before a timeout response or cancellation returns.
        # Replacement startup does not extend the failed job's HTTP deadline.
        termination = asyncio.create_task(self._stop(worker))
        try:
            await asyncio.shield(termination)
        except asyncio.CancelledError:
            # A disconnect may arrive during timeout cleanup. Finish reaping
            # before propagating that cancellation; still recover the capacity.
            await asyncio.shield(termination)
            raise
        finally:
            if (
                not self.closing and termination.done() and not termination.cancelled()
                and termination.exception() is None
            ):
                task = asyncio.create_task(self._recover(worker))
                self.recoveries.add(task)
                task.add_done_callback(self._recovered)

    def _recovered(self, task):
        self.recoveries.discard(task)
        if not task.cancelled() and task.exception() is not None:
            _internal_error(task.exception())

    async def _recover(self, worker):
        delay = 1
        while not self.closing:
            try:
                await self._start(worker)
                return
            except Exception as exc:
                await self._stop(worker)
                _internal_error(exc)
                await asyncio.sleep(delay)
                delay = min(delay * 2, 30)

    async def close(self):
        self.closing = True
        pending = [task for task in self.tasks if task is not asyncio.current_task()]
        for task in pending:
            task.cancel()
        if pending:
            await asyncio.gather(*pending, return_exceptions=True)
        recovering = list(self.recoveries)
        for task in recovering:
            task.cancel()
        if recovering:
            await asyncio.gather(*recovering, return_exceptions=True)
        await asyncio.gather(*(self._stop(worker) for worker in self.workers))

    async def execute(
        self, *, method, path, query=b"", body=b"", request_id=None,
        content_type="application/json",
    ):
        if (method, path) not in CHEMISTRY_ROUTES:
            raise ExecutionFailure("Operation is not a chemistry job", 400)
        if len(body) > self.config.max_request_bytes:
            raise ExecutionFailure("Request exceeds the production byte limit", 413)
        if len(query) > self.config.max_query_bytes:
            raise ExecutionFailure("Query exceeds the production byte limit", 413)
        if len(content_type) > 256:
            raise ExecutionFailure("Content-Type exceeds its limit", 413)
        worker = next((item for item in self.workers
                       if item.ready and not item.busy), None)
        if self.closing or worker is None:
            raise ExecutionFailure("Chemistry capacity is busy; retry shortly")
        worker.busy = True
        retired = False
        task = asyncio.current_task()
        self.tasks.add(task)
        message = {
            "method": method, "path": path,
            "query": base64.b64encode(query).decode("ascii"),
            "body": base64.b64encode(body).decode("ascii"),
            "request_id": request_id or uuid.uuid4().hex,
            "content_type": content_type,
        }

        async def exchange():
            data = _encode(message)
            await asyncio.get_running_loop().sock_sendall(
                worker.socket, struct.pack("!I", len(data)) + data
            )
            return await _read_async(
                worker.socket, self.config.max_response_bytes * 2 + 16384
            )

        async def retire_once():
            nonlocal retired
            # Cancellation can arrive while a result-triggered retirement is
            # awaiting process exit. That retirement already owns recovery.
            if not retired:
                retired = True
                await self._retire(worker)

        try:
            try:
                result = await asyncio.wait_for(exchange(), self.config.timeout_seconds)
                output = base64.b64decode(result["body"], validate=True)
                if len(output) > self.config.max_response_bytes:
                    raise ValueError("Worker response exceeds its limit")
                for diagnostic in result.get("diagnostics", [])[:8]:
                    _log("internal_error", **diagnostic)
                worker.jobs += 1
                if worker.jobs >= self.config.max_jobs or result["status"] >= 500:
                    await retire_once()
                return result["status"], result["headers"], output
            except asyncio.TimeoutError:
                await retire_once()
                _log("job_timeout", deadline_seconds=self.config.timeout_seconds)
                raise ExecutionFailure("Chemistry exceeded its time limit", 504)
            except asyncio.CancelledError:
                await retire_once()
                _log("job_cancelled")
                raise
            except ExecutionFailure:
                raise
            except Exception as exc:
                await retire_once()
                _internal_error(exc)
                raise ExecutionFailure("Chemistry worker restarted; retry the request")
        finally:
            worker.busy = False
            self.tasks.discard(task)


class ExecutionMiddleware:
    def __init__(self, app, state, config):
        self.app, self.state, self.config = app, state, config

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        route = scope["method"], scope["path"]
        registration = route == ("POST", "/register_monomer")
        if route not in CHEMISTRY_ROUTES and not registration:
            await self.app(scope, receive, send)
            return
        try:
            if registration:
                from .security import require_registration

                # FastAPI otherwise parses the entire model before entering the
                # endpoint guard. Reject unauthorized uploads without reading.
                require_registration(Request(scope))
            else:
                executor = getattr(self.state, "executor", None)
                if executor is None or not self.state.readiness.get("ready"):
                    raise ExecutionFailure("Chemistry is not ready")
                if not executor.available:
                    raise ExecutionFailure("Chemistry capacity is busy; retry shortly")
            body = bytearray()

            async def collect_body():
                while True:
                    event = await receive()
                    if event["type"] == "http.disconnect":
                        scope["cabiln.disconnected"] = True
                        raise asyncio.CancelledError()
                    chunk = event.get("body", b"")
                    if len(body) + len(chunk) > self.config.max_request_bytes:
                        raise ExecutionFailure("Request exceeds the byte limit", 413)
                    body.extend(chunk)
                    if not event.get("more_body"):
                        return

            try:
                await asyncio.wait_for(collect_body(), 10)
            except asyncio.TimeoutError:
                raise ExecutionFailure("Request upload exceeded its time limit", 408)

            if registration:
                replayed = False

                async def replay():
                    nonlocal replayed
                    if not replayed:
                        replayed = True
                        return {"type": "http.request", "body": bytes(body),
                                "more_body": False}
                    return await receive()

                # Authentication and byte admission finish before ordinary
                # model validation; the atomic administrator write stays local.
                await self.app(scope, replay, send)
                return

            async def disconnected():
                while (await receive())["type"] != "http.disconnect":
                    pass
                scope["cabiln.disconnected"] = True

            content_type = next(
                (value.decode("latin1") for key, value in scope.get("headers", [])
                 if key.lower() == b"content-type"),
                "application/json",
            )

            job = asyncio.create_task(executor.execute(
                method=scope["method"], path=scope["path"],
                query=scope.get("query_string", b""), body=bytes(body),
                request_id=_request_id.get(),
                content_type=content_type,
            ))
            gone = asyncio.create_task(disconnected())
            try:
                done, _ = await asyncio.wait({job, gone},
                                             return_when=asyncio.FIRST_COMPLETED)
                if gone in done:
                    job.cancel()
                    await asyncio.gather(job, return_exceptions=True)
                    return
                status, headers, output = await job
                response = Response(output, status_code=status)
                response.raw_headers = [
                    (key.encode("latin1"), value.encode("latin1"))
                    for key, value in headers
                    if key.lower() not in ("transfer-encoding", "connection")
                ]
            finally:
                gone.cancel()
                if not job.done():
                    job.cancel()
                await asyncio.gather(gone, job, return_exceptions=True)
        except ExecutionFailure as exc:
            response = JSONResponse(
                {"error": str(exc), "request_id": _request_id.get()},
                status_code=exc.status,
                headers={"Retry-After": "1"} if exc.status == 503 else None,
            )
        except HTTPException as exc:
            response = error_response(exc)
        await response(scope, receive, send)


class RequestLoggingMiddleware:
    def __init__(self, app, version):
        self.app, self.version = app, version

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        identity = uuid.uuid4().hex
        token = _request_id.set(identity)
        started, status = time.monotonic(), 500
        response_started = False

        async def record(event):
            nonlocal status, response_started
            if event["type"] == "http.response.start":
                response_started = True
                status = event["status"]
                event["headers"] = [*event.get("headers", []),
                                    (b"x-request-id", identity.encode())]
            await send(event)

        try:
            try:
                await self.app(scope, receive, record)
            except Exception as exc:
                if response_started:
                    _internal_error(exc)
                else:
                    await error_response(exc)(scope, receive, record)
            except asyncio.CancelledError:
                if not response_started:
                    status = 499
                raise
        finally:
            route = scope.get("route")
            path = scope.get("path")
            label = getattr(route, "path", None)
            if (scope["method"], path) in CHEMISTRY_ROUTES:
                label = path
            if scope.get("cabiln.disconnected") and not response_started:
                status = 499
            _log("http_request", route=label or "unmatched", status=status,
                 method=scope["method"], release=self.version,
                 duration_ms=round((time.monotonic() - started) * 1000, 2))
            _request_id.reset(token)
