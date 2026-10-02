# Development tools

The application uses the installed package entry points. These are the maintained
ways to build, validate, register and display peptides:

| Entry point | Purpose |
| --- | --- |
| `cabiln` | Start the GUI and HTTP application |
| `run_pyPept` | Build structures, depictions and conformers |
| `pyPept-BILN-validate` | Validate notation |
| `pyPept-monomer-add` | Activate and register monomers |
| `pyPept.interfaces.monomer_pipeline` | Import HELM libraries, derive CSV records and build SDF libraries |
| `pyPept.smiles.convert_smiles` | Import a molecular structure with residue ownership and recognition status |
| `pyPept.inputs` | Interpret formats and verify notation conversion for application callers |

`live_renderer.py` remains a compatibility launcher. Its historical imports and
text-only crosslink helper remain available. New chemistry code should import
the package directly; dedicated tests cover the launcher separately.

The maintained validation commands are:

```sh
python -m pytest -m "not distribution and not fuzz"
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
python tools/check_examples.py
python -m pyPept.library_quality --output /tmp/library-quality-candidate.json
python tools/benchmarks/decomposition.py --permutations 1 --output /tmp/cabiln-benchmark
```

The [browser suite](../tests/browser/README.md) checks GUI interactions against an
isolated library. The [decomposition benchmark](../docs/decomposition-benchmark.md)
uses independent reference structures, partitions and attachment expectations.
Its external comparator and recorded results live in `benchmarks/`.

The [library compatibility audit](../docs/library-quality.md) uses the package's
current restoration and activation rules. It records exact per-definition
results, metadata differences and curated identity concerns. Review the candidate
report before changing the packaged baseline.

Two additional diagnostics answer different questions:

| Command | Scope |
| --- | --- |
| `python tools/audit_library.py` | CSV structure duplicates, name/CIP heuristics and CSV/SDF template consistency |
| `python tools/find_duplicates.py` | SDF template/core duplicates and case-insensitive symbol collisions |

Both print reports without changing the library. Duplicate structures can be
intentional. Review the name/CIP heuristics and chemistry before merging
definitions or changing stereochemistry.

`release_manifest.py` records the installed wheel and chemistry bindings during
the Docker build. `release_smoke.py` checks a running release image. See the
[deployment instructions](../docs/deployment.md) for their use.

## Generated tests

Generated tests use Hypothesis for notation, chemistry and library changes, and
fast-check for UI histories and browser interactions. From the repository root,
after installing `[dev,web]`, run the bounded campaign:

```bash
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/python python -m pytest -m fuzz \
  --hypothesis-show-statistics --timeout=300 --timeout-method=thread
npm ci --prefix tests/browser
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/frontend npm run fuzz:frontend --prefix tests/browser
(cd tests/browser && npx playwright install chromium)
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/browser npm run fuzz:browser --prefix tests/browser
```

Set `CABILN_FUZZ_PROFILE=deep` for longer campaigns; use `--timeout=1200` for the
Python run. CI runs bounded campaigns on pushes and pull requests, and the deep
profile nightly. Browser servers and library changes use temporary local copies.
These campaigns sample supported feature combinations; they do not exhaust the
possible peptides or editing histories.

`CABILN_FUZZ_CASE=browser-tutorial-lifecycle` selects generated Build and Swap
lesson histories. It varies preview use, viewport width, guide and panel closures,
Back/Next and Undo/Redo, checking the document against the expected edit after
each action.

Failures retain reduced inputs, seeds, dependency versions and replay information.
Python artifacts include source/library hashes, Hypothesis's example database
and reproduction decorator. Replay with the recorded versions and
`--hypothesis-seed=<seed>`, or apply the reported reproduction decorator to the
owning test. Select the recorded fast-check property with `CABILN_FUZZ_CASE`, then
use `CABILN_FUZZ_SEED`, `CABILN_FUZZ_PATH` and, for command histories,
`CABILN_FUZZ_REPLAY_PATH` to replay its failure. `CABILN_FUZZ_SCHEDULE` accepts the
saved task order for exact replay of a frontend scheduling failure.
Observation counters include shrinking and repeated examples; they do not count
unique peptides. A watchdog termination preserves
the latest eight Python observations in `active.json`, alongside the captured CI log.

## Historical scripts

Dataset repair recipes, monomer batches and scratch probes are preserved in
[commit c67df6f](https://github.com/anagnorisis2peripeteia/pyPept/tree/c67df6fdc57e875cf6001ad9e5178a22188eb7d1/tools).
To inspect an original recipe:

```sh
git show c67df6f:tools/add_monomers_batch6.py
```

The former `validate_monomers.py` and `full_library_roundtrip.py` diagnostics are
superseded by the versioned library audit above. The README check is now
`check_examples.py`.

The root `bench_cyclicpepedia.py` survey depended on removed normalization
helpers. Its source remains in that commit; `cyclicpepedia_structure.xlsx` and
recorded benchmark results remain available. The current decomposition benchmark
uses independently annotated cases and has a different scope.
