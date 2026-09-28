# Browser acceptance

These tests drive the local application through its visible controls and actual
HTTP chemistry. Every server uses a temporary monomer library. Registration
never changes packaged data. Tests close their servers and remove the copies.

Install the project's Python `[dev,web]` dependencies, then run from this directory:

```sh
npm ci
npx playwright install chromium
npm test
```

The suite covers building chains, caps, branches, disulfides, cycles and
scaffolds; panels; tile families; chip/SVG selection; insertion; nested atom
highlights; input formats; verification/uploads; registration; unknown-region
editing; canvas controls; exports; and delayed/failed responses.

`CABILN_PYTHON` selects an interpreter. Otherwise the repository's `.venv` is used
when available, falling back to `python3`. `CABILN_APP_ROOT` selects another
checkout for baseline comparisons. Explicit `PYTHONPATH` keeps that source
selection independent of editable installations.

`CABILN_BROWSER_CHANNEL=chrome` uses installed Chrome. Without it, Playwright's
Chromium is used. The viewport is fixed at 1440 × 1000, with one worker.
`CABILN_BROWSER_ARTIFACTS` selects the report directory; the default is
`cabiln-browser-acceptance` beneath the system temporary directory.

Captures include source text, occurrence/group tabs, selection outlines, panel
bounds, enabled controls, atom highlights, screenshots and downloaded structures.
Screenshots settle animations for capture only; live interaction timing is
unchanged. Failure traces and server logs accompany the JSON result report.
Generated artifacts are not committed.

For a refactor, run the same journeys against a frozen baseline and the current
source with separate artifact directories. Compare captured state and geometry,
then screenshots. A new correctness fix gets its own failing baseline test; a
known defect must not become a compatibility expectation.
