"""Exercise an isolated running release image; emits no submitted structures."""

import argparse
import json
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8000")
    args = parser.parse_args()

    def request(path, body=None, *, token=None, status=200):
        data = json.dumps(body).encode() if body is not None else None
        req = Request(args.url + path, data=data)
        if data is not None:
            req.add_header("content-type", "application/json")
        if token is not None:
            req.add_header("X-Cabiln-Library", token)
        try:
            response = urlopen(req, timeout=35)
        except HTTPError as error:
            if error.code != status:
                raise
            response = error
        with response:
            assert response.status == status, (path, response.status)
            assert response.headers["x-request-id"]
            return json.load(response)

    deadline = time.monotonic() + 60
    while True:
        try:
            ready = request("/ready")
            break
        except (URLError, TimeoutError, ConnectionError):
            if time.monotonic() >= deadline:
                raise
            time.sleep(0.2)
    assert request("/capabilities") == {"registration": False, "session_registration": True}
    source = "K.[G(4,2)]-A"
    rendered = request("/render", {"cabiln": source})
    assert "<svg" in rendered["svg"] and len(rendered["residues"]) == 3
    prepared = request(
        "/prepare_project",
        {
            "project": {
                "format": "cabiln-project",
                "version": 1,
                "document": {"text": source, "notation": "cabiln"},
            }
        },
    )
    assert request("/validate_project", prepared)["valid"]
    assert request("/monomers")
    preview = request("/preview_monomer", {"smiles": "NCC(=O)O"})
    entry = {key: preview[key] for key in ("chuckles", "chem_types", "leaving", "activation_policy")}
    entry.update(abbr="ReleaseGly", name="Temporary release check")
    snapshot = request("/session_library", {"monomers": [entry]})
    rendered = request("/render", {"cabiln": "ReleaseGly"}, token=snapshot["token"])
    assert [item["abbr"] for item in rendered["residues"]] == ["ReleaseGly"]
    request("/render", {"cabiln": "ReleaseGly"}, status=400)
    request("/register_monomer", entry, status=403)
    print(json.dumps({"ready": ready, "checks": 8, "status": "passed"}))


if __name__ == "__main__":
    main()
