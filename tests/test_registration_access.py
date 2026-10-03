"""Public onboarding does not grant permission to change the installed library."""

import base64

import pytest
from fastapi.testclient import TestClient

from pyPept.web.app import create_app

TOKEN = "test-only-administrative-token-123456789"


def authorization(password=TOKEN):
    encoded = base64.b64encode(f"admin:{password}".encode()).decode()
    return {"authorization": f"Basic {encoded}"}


@pytest.fixture
def protected(monkeypatch):
    monkeypatch.setenv("CABILN_REGISTRATION_TOKEN", TOKEN)
    with TestClient(create_app(allow_registration=True)) as client:
        yield client


def test_public_onboarding_keeps_installed_writes_authenticated(protected):
    for headers in ({}, authorization("wrong"), {"authorization": "Basic @@@"}):
        page = protected.get("/register", headers=headers)
        assert page.status_code == 200
        assert TOKEN not in page.text
        response = protected.post("/register_monomer", headers=headers, json={
            "abbr": "TestAuth", "name": "test", "chuckles": "invalid",
            "chem_types": {}, "leaving": {},
        })
        assert response.status_code == 401
        assert "Basic" in response.headers["www-authenticate"]
        assert TOKEN not in response.text
    assert protected.get("/register", headers=authorization()).status_code == 200


def test_remote_address_cannot_use_local_registration(monkeypatch):
    monkeypatch.delenv("CABILN_REGISTRATION_TOKEN", raising=False)
    with TestClient(
        create_app(allow_registration=True), client=("198.51.100.42", 1234)
    ) as client:
        assert client.get("/register").status_code == 200
        assert client.post("/register_monomer", json={
            "abbr": "TestAuth", "name": "test", "chuckles": "invalid",
            "chem_types": {}, "leaving": {},
        }).status_code == 403
    with TestClient(
        create_app(allow_registration=True), client=("127.0.0.1", 1234)
    ) as client:
        assert client.get("/register").status_code == 200


@pytest.mark.parametrize(
    "headers",
    [
        {"origin": "https://unrelated.example"},
        {"sec-fetch-site": "cross-site"},
        {"origin": "null"},
        {"origin": "https://[malformed"},
    ],
)
def test_cross_origin_write_rejected_before_chemistry(protected, headers):
    payload = {
        "abbr": "TestAuth",
        "name": "test",
        "chuckles": "invalid",
        "chem_types": {},
        "leaving": {},
    }
    response = protected.post(
        "/register_monomer", json=payload, headers={**authorization(), **headers}
    )
    assert response.status_code == 403


def test_write_requires_credentials_independently_of_page(protected):
    protected.get("/register", headers=authorization())
    response = protected.post(
        "/register_monomer",
        json={
            "abbr": "TestAuth",
            "name": "test",
            "chuckles": "invalid",
            "chem_types": {},
            "leaving": {},
        },
    )
    assert response.status_code == 401


def test_public_registration_configuration_fails_closed(monkeypatch, tmp_path):
    monkeypatch.setenv("CABILN_ENV", "production")
    monkeypatch.delenv("CABILN_REGISTRATION_TOKEN", raising=False)
    with pytest.raises(ValueError, match="administrative token"):
        create_app(allow_registration=True, execution_mode="local")
    monkeypatch.setenv("CABILN_REGISTRATION_TOKEN", TOKEN)
    monkeypatch.delenv("CABILN_MONOMER_LIBRARY", raising=False)
    monkeypatch.delenv("CABILN_LIBRARY_BACKUP_DIR", raising=False)
    with pytest.raises(ValueError, match="external library"):
        create_app(allow_registration=True, execution_mode="local")
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(tmp_path / "monomers.sdf"))
    monkeypatch.setenv("CABILN_LIBRARY_BACKUP_DIR", str(tmp_path / "backups"))
    assert create_app(allow_registration=True, execution_mode="local")


def test_weak_administrative_token_is_not_accepted(monkeypatch):
    monkeypatch.setenv("CABILN_REGISTRATION_TOKEN", "short")
    with pytest.raises(ValueError, match="32"):
        create_app(allow_registration=True)
