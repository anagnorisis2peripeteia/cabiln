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
intentional; name/CIP heuristics need chemical review. They do not authorize
merging definitions or changing stereochemistry.

`release_manifest.py` records the installed wheel and chemistry bindings during
the Docker build. `release_smoke.py` checks a running release image. See the
[deployment instructions](../docs/deployment.md) for their use.

## Historical scripts

Dataset repair recipes, monomer batches and scratch probes are preserved in
[commit 8721900](https://github.com/anagnorisis2peripeteia/pyPept/tree/872190083427bef49f2e8d94af827eea74080775/tools).
To inspect an original recipe:

```sh
git show 8721900:tools/add_monomers_batch6.py
```

The former `validate_monomers.py` and `full_library_roundtrip.py` diagnostics are
superseded by the versioned library audit above. The README check is now
`check_examples.py`. Past review reports retain their original paths and counts.

The root `bench_cyclicpepedia.py` survey depended on removed normalization
helpers. Its source remains in that commit; `cyclicpepedia_structure.xlsx` and
recorded benchmark results remain available. The current decomposition benchmark
uses independently annotated cases and has a different scope.
