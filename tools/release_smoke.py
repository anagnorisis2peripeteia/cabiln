"""Exercise an isolated running release image; emits no submitted structures."""

import argparse
import json
import time
from urllib.error import URLError
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8000")
    args = parser.parse_args()

    def request(path, body=None):
        data = json.dumps(body).encode() if body is not None else None
        req = Request(args.url + path, data=data)
        if data is not None:
            req.add_header("content-type", "application/json")
        with urlopen(req, timeout=35) as response:
            assert response.status == 200
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
    assert request("/capabilities") == {"registration": False}
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
    print(json.dumps({"ready": ready, "checks": 5, "status": "passed"}))


if __name__ == "__main__":
    main()
