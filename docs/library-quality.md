# Bundled monomer compatibility audit

`src/pyPept/data/library-quality.json` records the current software behavior for
each of the 1,128 bundled monomer definitions. It replaces the former aggregate
activation-success threshold with an exact per-record regression baseline.
It is not a certificate of scientific identity, stereochemical completeness,
reaction yield, or complete chemical validation.

## What is recorded

The audit starts from the bundled `monomers.sdf`, independently of any custom
library selected in the environment. For every record it stores:

- A definition hash covering canonical isomeric template SMILES, named identity,
  monomer type, numbered leaving groups, and declared site chemistry. Coordinates
  and atom traversal are excluded. The complete SDF also has a file hash.
- Expected and observed activation templates and leaving groups. Restoration
  uses the same shared leaving-group implementation as molecule assembly; the
  restored standalone molecule is then passed to automatic activation.
- Explicit issues for legacy R3 sidechain numbering, absent leaving metadata,
  differences between declared and effective site chemistry, unspecified
  tetrahedral stereochemistry, and activation exceptions.

Legacy R3 sidechain records are compared after shifting their sidechain slot
numbers by one in a temporary copy. Their original slots and structures remain
unchanged. This comparison accommodates the historical numbering convention;
it does not migrate a library or user document.

The current activation baseline is 989 exact matches, 64 matches under that
legacy comparison, 73 mismatches, and two activation failures (`Hsl`, `TATA`).
All 75 activation exceptions identify the particular record and discrepancy.
Issue counts overlap: 66 records use the legacy convention, four have sparse
leaving metadata, 47 contain unspecified tetrahedral centers, and 270 have a
declared/effective chemistry difference. A difference is a review fact, not
proof that the stored chemical structure is incorrect. No configuration is
invented for an unspecified stereocenter.

## Runtime meaning

`pyPept.library_quality.quality_for_monomer(molecule)` performs a read-only
fingerprint lookup, without running activation or modifying any library. The
palette exposes that result as `quality`:

| Status | Meaning |
| --- | --- |
| `no_known_exception` | This exact definition matches the baseline and has none of its recorded issues. |
| `review_required` | This exact definition matches the baseline and has explicit issues to inspect. |
| `unreviewed` | This custom, renamed, changed, or otherwise unmatched definition has no applicable baseline. |

An unreviewed definition remains usable. Its empty `issues` array must not be
read as a clean bill of health. A changed fingerprint cannot inherit the known
status of another record with the same label. Unsupported enhanced stereo also
cannot inherit a fingerprint produced by plain-SMILES normalization.

`audit_version` names the fingerprint/audit convention. The manifest records its
generating RDKit version for provenance; its definition hashes and all per-record
facts have also been compared on the supported Python 3.9 and 3.13 environments
with their installed RDKit versions. A future canonical-SMILES change that alters
a hash fails the audit regression and makes that runtime definition unreviewed
until checked. No cross-version scientific identity guarantee is implied.

## Reproducing and reviewing changes

Run the full audit into a candidate file before changing the packaged baseline:

```sh
python -m pyPept.library_quality --output /tmp/library-quality-candidate.json
python -m pytest tests/test_library_quality.py \
  tests/test_bond_validation_and_assembly.py::TestLibraryRoundTrip
```

Generation is deterministic for a fixed SDF, code, and RDKit version. The test
compares every record, all issue facts, activation results, and the SDF hash. A
new failure cannot be hidden by a high aggregate success percentage. Review the
individual expected/observed templates and leaving groups before accepting any
manifest change; regenerating the file by itself is not chemical review.

The baseline was established with the shared restoration implementation and
explicit review of the exception families. Independently written structure
controls cover legacy N-methylcysteine and its disulfide, sparse-metadata Mpa
and its disulfide, unspecified aMeLeu stereochemistry, the standalone Hsl and
TATA products, and lysine sidechain acylation. These controls do not prove the
identity of every imported record. Source attribution and scientific curation
of individual exceptions remain separate work.

## Input normalization boundary

Reference and MOL rendering reject enhanced AND/OR stereo before any conversion
to plain SMILES. The same shared support guard is used by direct conversion.
Ordinary and absolute tetrahedral stereo remain supported, with exact-structure
controls through normalization, recognition, and verification. This prevents a
normalized reference from appearing exact after its original relative/mixture
stereo meaning was discarded.

Existing library symbols use the selection schema, which admits punctuation
already present in the bundled library; new registration names retain their
stricter naming rule. Selected slot numbers must be positive and resolve to
real attachment sites. They have no arbitrary R64 cutoff; request and execution
budgets are enforced separately from chemical slot numbering.
