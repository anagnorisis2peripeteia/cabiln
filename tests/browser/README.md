# Browser acceptance

These tests drive the local application through its visible controls and actual
HTTP chemistry. Every server uses a temporary monomer library. Registration
never changes packaged data. Tests close their servers and remove the copies.

Install the project's Python `[dev,web]` dependencies and Node 20 or newer,
then run from this directory. The command below exercises the process executor
used in production; omit `CABILN_EXECUTION=process` for local thread-pool mode:

```sh
npm ci
npx playwright install chromium firefox webkit
CABILN_EXECUTION=process npm test
CABILN_EXECUTION=process CABILN_BROWSER_ENGINE=firefox npm test
CABILN_EXECUTION=process CABILN_BROWSER_ENGINE=webkit npm test
```

The suite covers building chains, caps, branches, disulfides, cycles and
scaffolds; connection/swap previews, site mappings and attachment highlights; panels; tile
families; chip/SVG selection; insertion; nested atom highlights; input formats;
verification/uploads; registration; project Save/Open; unknown-region editing; canvas controls; exports; and delayed/failed responses.

The UX journeys also cover retained drawings and zoom, disabled stale selections,
keyboard building, attachment highlights, Undo/Redo, format drafts, recovery,
unavailable storage, preserved stereo warnings, library revision refreshes, and
building/verification at a 390 px viewport. Interaction journeys cover independent
Escape dismissal and focus return, keyboard examples, shared zoom/reset controls,
touch selection, interrupted dragging, in-place Retry, and visible warnings with
expandable details. Chromium additionally exercises real multi-finger pan/pinch
through CDP; that input API is unavailable in the other engines.
Registration checks empty input, failed previews, duplicate abbreviations and
keyboard submission. Draft clearing checks recovery before a new edit.
The accessibility journey runs axe-core on empty, library, selected residue,
Swap, Help, Examples, preview, Verify, light drawing, tutorial, narrow layout and
registration states. Reduced motion keeps contrast checks out of transitional
colours. Reports retain incomplete checks for manual review; a passing scan is
not an accessibility certification or a physical-device/screen-reader test.
Library journeys check actual matches when the selected site changes, filter
reset and search, and recovery from missing or stalled reaction data. Tutorial journeys
exercise the full build/preview/Undo route at desktop and mobile widths, failed
drawings, keyboard dismissal, restart and isolation from saved browser drafts.
They also cover reopening closed panels, browsing without search, changing the
chosen amino acid, and invalidating unfinished steps when selections are cleared.
Direct Connect, completion while the guide is closed, Back after Apply and Undo,
and Show control on a scrolled library have regression journeys. The Swap lesson
loads Retatrutide, reviews its three occupied sites, replaces K17 and undoes the
edit at desktop and mobile widths. Example loading covers failure and switching
lessons while a response is pending.
Selection guidance checks the exact owned atoms, visibility without hover, pointer
selection, Highlight off, dismissal, Back and source changes followed by Undo.
Delayed previews keep a waiting instruction; review text follows the completed
edit even after choosing another block in Library.
Drawing selection checks bond margins and label interiors, zoom, touch, stale
drawings, and rejection of background clicks and pan gestures.

For comparable browser timings, run the performance file alone against each
revision with the same Python and browser:

```sh
CABILN_BROWSER_ARTIFACTS=/tmp/cabiln-current-timing npx playwright test performance.spec.js
CABILN_APP_ROOT=/path/to/frozen/source CABILN_BROWSER_ARTIFACTS=/tmp/cabiln-baseline-timing npx playwright test performance.spec.js
```

Its `browser-timings.json` records `feedbackMs` (the first visible acknowledgement)
and `completedMs` (the finished drawing and chips), each after two animation
frames. It covers edits, connections, a 60-residue chain with caps and six nested
branches. Completion includes network and chemistry; assertion polling is excluded.
The first edit uses a fresh server only when this file runs alone. These are small
local samples, not field INP or production percentiles. The suite requires actual
feedback, with no machine-specific speed threshold.

`CABILN_PYTHON` selects an interpreter. Otherwise the repository's `.venv` is used
when available, falling back to `python3`. `CABILN_APP_ROOT` selects another
checkout for baseline comparisons. Explicit `PYTHONPATH` keeps that source
selection independent of editable installations. `CABILN_BROWSER_INSTALLED=1`
uses the installed wheel with Python isolation; release CI sets this and runs
with the production worker memory limit.

`CABILN_BROWSER_CHANNEL=chrome` uses installed Chrome. Without it, Playwright's
Chromium is used. `CABILN_BROWSER_ENGINE=firefox` or `webkit` selects another
engine. CI runs the suite on all three against the same release wheel. The
default viewport is 1440 × 1000, with one worker; responsive journeys override it.
On macOS, WebKit follows the system's keyboard navigation preference. Enable
Keyboard navigation before testing Tab between buttons, or set
`defaults write org.webkit.Playwright AppleKeyboardUIMode -int 2` for the test
browser only. Restore its previous value after testing; use `defaults delete
org.webkit.Playwright AppleKeyboardUIMode` if no value was previously set.
The skip-link journey uses Option+Tab on macOS WebKit to include links, matching
[Safari's keyboard behaviour](https://github.com/microsoft/playwright/issues/5609).
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
