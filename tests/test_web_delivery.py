"""Delivery optimizations preserve assets and molecular render payloads."""

import gzip

from fastapi.testclient import TestClient
import pytest

from pyPept.web.app import create_app


@pytest.mark.parametrize(
    "source,word",
    [
        ("A-G-Nope-G", "Nope"),
        ("K.[G(4,2).Missing(1,2)]-A", "Missing"),
        ("A-🧬-G", "🧬"),
    ],
)
def test_render_errors_keep_original_source_locations_and_repair_hints(source, word):
    with TestClient(create_app()) as client:
        response = client.post("/render", json={"cabiln": source})
    assert response.status_code == 400
    diagnostic = response.json()
    assert "Library" in diagnostic["hint"]
    span = diagnostic["source_span"]
    assert source[span["start"] : span["end"]] == word
    assert word in diagnostic["error"]


def test_error_in_rewritten_legacy_source_does_not_point_into_the_wrong_text():
    with TestClient(create_app()) as client:
        response = client.post("/render", json={"cabiln": "A-D.(4,1)-A-am%G-Nope-am"})
    assert response.status_code == 400
    diagnostic = response.json()
    assert "Nope" in diagnostic["error"]
    assert "Library" in diagnostic["hint"]
    assert "source_span" not in diagnostic


@pytest.mark.parametrize("path", ["/", "/static/builder.js", "/static/builder.css"])
def test_browser_assets_negotiate_compression_without_changing_content(path):
    with TestClient(create_app()) as client:
        plain = client.get(path, headers={"Accept-Encoding": "identity"})
        with client.stream("GET", path, headers={"Accept-Encoding": "gzip"}) as compressed:
            transferred = b"".join(compressed.iter_raw())
    assert plain.status_code == compressed.status_code == 200
    assert "content-encoding" not in plain.headers
    assert compressed.headers["content-encoding"] == "gzip"
    assert "Accept-Encoding" in compressed.headers["vary"]
    # File responses can use chunked transfer without a Content-Length header.
    assert gzip.decompress(transferred) == plain.content
    assert len(transferred) < len(plain.content) / 2


def test_compressed_render_keeps_svg_mol_and_atom_ownership():
    with TestClient(create_app()) as client:
        body = {"cabiln": "ac-K.[G(4,2).A(1,2)]-am", "width": 900, "height": 600}
        plain = client.post("/render", json=body, headers={"Accept-Encoding": "identity"})
        compressed = client.post("/render", json=body, headers={"Accept-Encoding": "gzip"})
        health = client.get("/health", headers={"Accept-Encoding": "gzip"})
    assert plain.status_code == compressed.status_code == 200
    assert compressed.headers["content-encoding"] == "gzip"
    assert compressed.json() == plain.json()
    assert int(compressed.headers["content-length"]) < len(plain.content) / 2
    assert "content-encoding" not in health.headers
