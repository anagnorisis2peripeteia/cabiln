"""Production request admission, actual process lifetimes and safe diagnostics."""

import asyncio
import base64
from dataclasses import replace
import importlib
import json
import logging
import os
from pathlib import Path
import sys
import time

from fastapi.testclient import TestClient
import pytest

from pyPept.web import cache
from pyPept.web.app import create_app
from pyPept.web.execution import (
    ChemistryExecutor, ExecutionConfig, ExecutionFailure, error_response,
)


def test_startup_readiness_prepares_selected_library_once(monkeypatch):
    module = importlib.import_module("pyPept.web.app")
    original, calls = module.check_readiness, []
    encode, palettes = module.JSONResponse.render, []

    def checked(directory):
        calls.append(directory)
        return original(directory)

    def encoded(response, content):
        if isinstance(content, list) and content and "abbr" in content[0]:
            palettes.append(content)
        return encode(response, content)

    monkeypatch.setattr(module, "check_readiness", checked)
    monkeypatch.setattr(module.JSONResponse, "render", encoded)
    cache.clear_render_cache()
    with TestClient(create_app()) as client:
        # Readiness includes chemistry metadata and its JSON representation;
        # the first user must not prepare the whole palette on their request.
        assert len(palettes) == 1
        initial = client.get("/monomers")
        assert initial.status_code == 200
        assert initial.json() == palettes[0]
        for _ in range(3):
            assert client.get("/ready").status_code == 200
            assert client.get("/health").json() == {"status": "ok"}
            current = client.get("/monomers")
            assert current.content == initial.content
            assert current.headers["x-library-version"] == initial.headers["x-library-version"]
    assert len(calls) == 1
    assert len(palettes) == 1


@pytest.mark.parametrize("broken", ["missing", "invalid", "assets", "palette"])
def test_readiness_fails_without_usable_data_or_assets(monkeypatch, tmp_path, broken):
    if broken == "palette":
        module = importlib.import_module("pyPept.web.monomers")

        def fail():
            raise RuntimeError("private-library palette failure")

        monkeypatch.setattr(module, "palette_response", fail)
    elif broken == "assets":
        module = importlib.import_module("pyPept.web.app")
        monkeypatch.setattr(module, "STATIC_DIR", tmp_path)
    else:
        library = tmp_path / "private-library.sdf"
        if broken == "invalid":
            library.write_text("malformed record\n$$$$\n", encoding="utf-8")
        monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(library))
    with TestClient(create_app()) as client:
        ready = client.get("/ready")
        assert ready.status_code == 503
        assert ready.json()["status"] == "not_ready"
        assert "private-library" not in ready.text
        assert client.get("/health").json() == {"status": "ok"}


def test_readiness_requires_the_shipped_quality_manifest(monkeypatch):
    module = importlib.import_module("pyPept.web.readiness")
    monkeypatch.setattr(module, "quality_manifest", lambda: {"schema_version": 99})
    with TestClient(create_app()) as client:
        assert client.get("/ready").status_code == 503
        assert client.get("/health").status_code == 200


def test_release_is_deployment_wide_and_diagnostics_do_not_log_input(
    monkeypatch, caplog,
):
    monkeypatch.setenv("CABILN_RELEASE", "reviewed-release")
    app = create_app()

    @app.get("/runtime-failure")
    def fail():
        try:
            raise RuntimeError("SECRET-MOLECULE-C[C@@H](N)C(=O)O")
        except Exception as exc:
            return error_response(exc)

    with caplog.at_level(logging.INFO, logger="cabiln.runtime"):
        with TestClient(app) as client:
            response = client.get("/runtime-failure?smiles=PRIVATE-QUERY")
            assert client.get("/server_id").text == "reviewed-release"
        with TestClient(create_app()) as client:
            assert client.get("/server_id").text == "reviewed-release"
    assert response.status_code == 500
    assert response.json() == {
        "error": "An internal error occurred.",
        "request_id": response.headers["x-request-id"],
    }
    records = [record.message for record in caplog.records
               if record.name == "cabiln.runtime"]
    text = "\n".join(records)
    assert "SECRET-MOLECULE" not in text
    assert "PRIVATE-QUERY" not in text
    assert "RuntimeError" in text
    assert '"duration_ms":' in text
    assert response.headers["x-request-id"] in text


def test_uncaught_endpoint_failure_has_safe_request_id(caplog):
    app = create_app()

    @app.get("/uncaught-runtime-failure")
    async def fail():
        raise RuntimeError("UNSAFE-EXCEPTION-TEXT")

    with TestClient(app) as client:
        response = client.get("/uncaught-runtime-failure")
    assert response.status_code == 500
    assert response.json()["request_id"] == response.headers["x-request-id"]
    assert "UNSAFE-EXCEPTION-TEXT" not in caplog.text


def test_render_cache_accounts_bytes_overwrite_lru_and_oversize(monkeypatch):
    cache.clear_render_cache()
    value = {"svg": "x" * 1000, "residues": [1, 2, 3]}
    size = cache._retained_size(("a", value))
    monkeypatch.setattr(cache, "_RENDER_CACHE_BYTES", size * 2)
    try:
        cache._rc_put("a", value)
        cache._rc_put("b", value)
        assert cache._rc_get("a") == value
        cache._rc_put("c", value)
        assert cache._rc_get("b") is None
        assert cache._cache_bytes <= size * 2
        before = cache._cache_bytes
        cache._rc_put("a", value)
        assert cache._cache_bytes == before
        cache._rc_put("a", {"svg": "y" * (size * 3)})
        assert cache._rc_get("a") is None
        assert cache._cache_bytes == size
    finally:
        cache.clear_render_cache()


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


async def _wait_recovered(pool):
    async def ready():
        while not pool.ready:
            await asyncio.sleep(0.02)
    await asyncio.wait_for(ready(), 20)


async def _wait_busy(pool):
    async def busy():
        while not pool.workers[0].busy:
            await asyncio.sleep(0.001)
    await asyncio.wait_for(busy(), 5)


async def _render(pool, source="K.[G(4,2).[ac(1,2)]]-A"):
    return await pool.execute(
        method="POST", path="/render",
        body=json.dumps({"cabiln": source}).encode(),
    )


def test_real_process_timeout_overload_cancellation_recovery_and_shutdown():
    async def exercise():
        pool = ChemistryExecutor(ExecutionConfig(mode="process", timeout_seconds=20))
        pids = []
        try:
            await pool.start()
            pids.append(pool.workers[0].process.pid)
            status, _, body = await _render(pool)
            assert status == 200
            assert json.loads(body)["bracket_groups"]

            # A deadline applies to real chemistry. The long peptide cannot
            # finish within 1 ms, and its actual process must die before504.
            pool.config = replace(pool.config, timeout_seconds=0.001)
            active = asyncio.create_task(_render(pool, "-".join(["G"] * 200)))
            await asyncio.sleep(0)
            started = time.monotonic()
            with pytest.raises(ExecutionFailure) as busy:
                await _render(pool)
            assert busy.value.status == 503
            assert time.monotonic() - started < 0.5
            with pytest.raises(ExecutionFailure) as timed_out:
                await active
            assert timed_out.value.status == 504
            assert time.monotonic() - started < 3
            assert not _alive(pids[-1])
            pool.config = replace(pool.config, timeout_seconds=20)
            await _wait_recovered(pool)
            pids.append(pool.workers[0].process.pid)
            assert pids[-1] != pids[-2]

            active = asyncio.create_task(_render(pool, "-".join(["G"] * 200)))
            await asyncio.sleep(0.005)
            active.cancel()
            with pytest.raises(asyncio.CancelledError):
                await active
            assert not _alive(pids[-1])
            await _wait_recovered(pool)
            pids.append(pool.workers[0].process.pid)
            assert (await _render(pool))[0] == 200

            # A dead child is replaced when transport detects the lost worker.
            pool.workers[0].process.kill()
            with pytest.raises(ExecutionFailure) as crashed:
                await _render(pool)
            assert crashed.value.status == 503
            await _wait_recovered(pool)
            pids.append(pool.workers[0].process.pid)
            assert (await _render(pool))[0] == 200
        finally:
            await pool.close()
        assert all(not _alive(pid) for pid in pids)
        assert not pool.ready
        assert not pool.tasks and not pool.recoveries
    asyncio.run(exercise())


def test_cancellation_during_retirement_cannot_start_two_replacements(monkeypatch):
    async def exercise():
        module = importlib.import_module("pyPept.web.execution")
        pool = ChemistryExecutor(ExecutionConfig(mode="process"))
        pool.workers[0].ready = True
        stopping, release = asyncio.Event(), asyncio.Event()
        calls = {"stop": 0, "recover": 0}

        async def sent(*args):
            pass

        async def response(*args):
            return {"status": 500, "headers": [], "body": "e30="}

        async def stop(worker):
            worker.ready = False
            calls["stop"] += 1
            stopping.set()
            await release.wait()

        async def recover(worker):
            calls["recover"] += 1

        monkeypatch.setattr(asyncio.get_running_loop(), "sock_sendall", sent)
        monkeypatch.setattr(module, "_read_async", response)
        monkeypatch.setattr(pool, "_stop", stop)
        monkeypatch.setattr(pool, "_recover", recover)
        active = asyncio.create_task(_render(pool))
        await stopping.wait()
        active.cancel()
        release.set()
        with pytest.raises(asyncio.CancelledError):
            await active
        await asyncio.gather(*pool.recoveries)
        assert calls == {"stop": 1, "recover": 1}
        assert not pool.tasks
    asyncio.run(exercise())


def test_process_transport_bounds_routes_and_complete_http_validation():
    async def exercise():
        config = ExecutionConfig(mode="process", max_response_bytes=1024)
        pool = ChemistryExecutor(config)
        try:
            # Admission rejects these before allocating a worker or IPC message.
            with pytest.raises(ExecutionFailure) as forbidden:
                await pool.execute(method="POST", path="/register_monomer")
            assert forbidden.value.status == 400
            with pytest.raises(ExecutionFailure) as oversized:
                await pool.execute(
                    method="POST", path="/render",
                    body=b"x" * (config.max_request_bytes + 1),
                )
            assert oversized.value.status == 413
            await pool.start()
            status, _, body = await pool.execute(
                method="POST", path="/render", body=b'{"cabiln":"G","width":0}',
            )
            assert status == 422
            assert "width" in json.loads(body)["error"]
            status, _, _ = await _render(pool, "G")
            assert status == 413
            assert pool.ready
        finally:
            await pool.close()
    asyncio.run(exercise())


def test_process_app_health_stays_responsive_while_worker_is_busy(monkeypatch):
    monkeypatch.setenv("CABILN_EXECUTION", "process")
    monkeypatch.setenv("CABILN_WORKER_MEMORY_MB", "0")
    app = create_app()
    with TestClient(app) as client:
        assert client.get("/ready").status_code == 200
        # Occupied admission is deterministic; the real lifecycle test above
        # verifies active chemistry, termination and replacement independently.
        app.state.executor.workers[0].busy = True
        try:
            started = time.monotonic()
            assert client.get("/health").json() == {"status": "ok"}
            assert client.get("/ready").status_code == 200
            response = client.post("/render", json={"cabiln": "G"})
            assert response.status_code == 503
            assert response.headers["retry-after"] == "1"
            assert time.monotonic() - started < 1
        finally:
            app.state.executor.workers[0].busy = False
        rendered = client.post("/render", json={"cabiln": "G"})
        assert rendered.status_code == 200
        assert rendered.headers["x-request-id"]


def test_http_disconnect_and_shutdown_terminate_active_chemistry(monkeypatch):
    monkeypatch.setenv("CABILN_WORKER_MEMORY_MB", "0")

    async def invoke(app, path, payload=None, disconnect=None):
        events, sent = [], False
        scope = {
            "type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1",
            "method": "POST" if payload else "GET", "scheme": "http",
            "path": path, "raw_path": path.encode(), "query_string": b"",
            "root_path": "", "headers": [(b"content-type", b"application/json")],
            "client": ("127.0.0.1", 0), "server": ("test", 0),
        }

        async def receive():
            nonlocal sent
            if not sent:
                sent = True
                return {"type": "http.request", "body": json.dumps(payload).encode(),
                        "more_body": False}
            if disconnect is not None:
                await disconnect.wait()
                return {"type": "http.disconnect"}
            await asyncio.Future()

        async def send(event):
            events.append(event)

        await app(scope, receive, send)
        return events

    async def exercise():
        app = create_app(execution_mode="process")
        active = None
        async with app.router.lifespan_context(app):
            pool = app.state.executor
            assert pool.ready
            pid = pool.workers[0].process.pid
            disconnect = asyncio.Event()
            payload = {"cabiln": "-".join(["G"] * 200)}
            active = asyncio.create_task(invoke(app, "/render", payload, disconnect))
            await _wait_busy(pool)
            health = await asyncio.wait_for(invoke(app, "/health"), 1)
            assert health[0]["status"] == 200
            busy = await invoke(app, "/render", {"cabiln": "G"})
            assert busy[0]["status"] == 503
            disconnect.set()
            assert await asyncio.wait_for(active, 3) == []
            assert not _alive(pid)
            await _wait_recovered(pool)
            pid = pool.workers[0].process.pid
            active = asyncio.create_task(invoke(app, "/render", payload))
            await _wait_busy(pool)
        outcome = await asyncio.gather(active, return_exceptions=True)
        assert isinstance(outcome[0], asyncio.CancelledError)
        assert not _alive(pid)
        assert not pool.tasks and not pool.recoveries
    asyncio.run(exercise())


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux RLIMIT_AS")
def test_linux_worker_memory_limit_supports_real_chemistry_and_is_enforced():
    async def exercise():
        pool = ChemistryExecutor(ExecutionConfig(mode="process", memory_mb=1024))
        try:
            await pool.start()
            limits = Path(f"/proc/{pool.workers[0].process.pid}/limits").read_text()
            address_limit = next(line for line in limits.splitlines()
                                 if line.startswith("Max address space"))
            assert address_limit.split()[3:5] == [str(1024 ** 3)] * 2
            assert (await _render(pool))[0] == 200
        finally:
            await pool.close()
        insufficient = ChemistryExecutor(ExecutionConfig(mode="process", memory_mb=1))
        with pytest.raises(Exception):
            await insufficient.start()
        assert not insufficient.ready
        assert all(worker.process is None for worker in insufficient.workers)
    asyncio.run(exercise())


def test_production_requires_execution_and_memory_limits(monkeypatch):
    monkeypatch.setenv("CABILN_ENV", "production")
    monkeypatch.setenv("CABILN_EXECUTION", "local")
    with pytest.raises(ValueError, match="requires CABILN_EXECUTION"):
        create_app()
    monkeypatch.setenv("CABILN_EXECUTION", "process")
    monkeypatch.setenv("CABILN_WORKER_MEMORY_MB", "0")
    with pytest.raises(ValueError, match="positive worker memory"):
        create_app()


def test_cli_refuses_unauthenticated_registration_on_public_bind(monkeypatch):
    module = importlib.import_module("pyPept.web.app")
    monkeypatch.delenv("CABILN_REGISTRATION_TOKEN", raising=False)
    monkeypatch.setattr(module, "app", create_app(allow_registration=False))
    monkeypatch.setattr(sys, "argv", ["cabiln", "--host", "0.0.0.0",
                                     "--enable-registration"])
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2


def test_registration_auth_precedes_body_read_and_authorized_body_is_bounded(
    monkeypatch,
):
    secret = "administrative-test-secret-32-characters"
    monkeypatch.setenv("CABILN_REGISTRATION_TOKEN", secret)
    monkeypatch.setenv("CABILN_MAX_REQUEST_BYTES", "1024")
    app = create_app(allow_registration=True, execution_mode="process")

    async def invoke(chunks, *, authenticated):
        headers = [(b"content-type", b"application/json"), (b"host", b"test")]
        if authenticated:
            credentials = base64.b64encode(f"admin:{secret}".encode())
            headers.append((b"authorization", b"Basic " + credentials))
        scope = {
            "type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1",
            "method": "POST", "scheme": "http", "path": "/register_monomer",
            "raw_path": b"/register_monomer", "query_string": b"",
            "root_path": "", "headers": headers,
            "client": ("203.0.113.1", 0), "server": ("test", 0),
        }
        events, consumed = [], []

        async def receive():
            if not authenticated:
                pytest.fail("Unauthorized registration body must not be read")
            consumed.append(chunks.pop(0))
            return {"type": "http.request", "body": consumed[-1],
                    "more_body": bool(chunks)}

        async def send(event):
            events.append(event)

        await app(scope, receive, send)
        return events[0]["status"]

    assert asyncio.run(invoke([], authenticated=False)) == 401
    assert asyncio.run(invoke([b"x" * 800, b"y" * 800], authenticated=True)) == 413
    assert asyncio.run(invoke([b"{}"], authenticated=True)) == 422
