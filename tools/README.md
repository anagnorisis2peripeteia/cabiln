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
python -m pytest -m "not distribution"
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
python tools/benchmarks/decomposition.py --permutations 1 --output /tmp/cabiln-benchmark.json
```

The [browser suite](../tests/browser/README.md) checks GUI interactions against an
isolated library. The [decomposition benchmark](../docs/decomposition-benchmark.md)
uses independent reference structures, partitions and attachment expectations.
Its external comparator and recorded results live in `benchmarks/`.

The remaining scripts are retained development history and library-maintenance
utilities. They are not imported by the application or part of its startup:

- `add_monomers_batch*.py`, `add_new_monomers.py` and `add_dye_monomers.py`
  record historical library additions.
- `fix_*.py`, `patch_monomers_*.py`, `update_sdf*.py`, `populate_chem_types.py`
  and `dedup_library.py` record data migrations and repairs.
- `audit_library.py`, `validate_monomers.py`, `full_library_roundtrip.py`,
  `check_*.py`, `find_duplicates.py`, `verify_chains.py`, `review_conflicts.py`,
  `investigate_fewer_slots.py` and related analysis scripts inspect past datasets.
- Files beginning with `_` and the standalone `test_*.py` probes are historical
  experiments. Current regression checks live under `tests/`.

These callable paths and repair capabilities remain in place. Several scripts
assume an old working directory or write directly to the bundled library; inspect
their paths and use a disposable library copy when reproducing a past migration.
Their name-specific repairs do not define runtime attachment or recognition rules.
