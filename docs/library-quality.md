# Bundled monomer compatibility audit

`src/pyPept/data/library-quality.json` records the current software behavior for
each of the 1,128 bundled monomer definitions. It replaces the former aggregate
activation-success threshold with an exact per-record regression baseline.
It is not a certificate of scientific identity, stereochemical completeness,
reaction yield, or complete chemical validation.

## What is recorded

The audit starts from the bundled `monomers.sdf`, independently of any custom
library selected in the environment. For every record it stores:

- A definition hash covering canonical isomeric template SMILES, symbol, full name,
  monomer type, numbered leaving groups, and declared site chemistry. Coordinates
  and atom traversal are excluded. The complete SDF also has a file hash.
- Expected and observed activation templates and leaving groups. Restoration
  uses the same shared leaving-group implementation as molecule assembly; the
  restored standalone molecule is then passed to automatic activation.
- Issues with an explicit `info` or `warning` severity for legacy numbering,
  absent or orphaned metadata, unresolved site chemistry, unspecified
  tetrahedral stereochemistry, activation exceptions, and reviewed identity
  conflicts. Each warning states the affected fact.

Legacy R3 sidechain records are compared after shifting their sidechain slot
numbers by one in a temporary copy. Their original slots and structures remain
unchanged. This comparison accommodates the historical numbering convention;
it does not migrate a library or user document.

The current activation baseline is 951 exact matches, 64 matches under that
legacy comparison, 101 mismatches, ten definitions requiring an orientation
choice, and two activation failures (`Hsl`, `TATA`).
The audit records every expected and observed template. It also checks whether
reactivation preserves the complete standalone isomeric structure. Preserved
structures with different selected sites or numbering produce informational
notes. Legacy numbering alone also produces a note.
An orientation choice is informational: the stored definition already states
its intended sites. Its audit record lists the possible reactivation templates.
Canonical R4+ numbering and detection of secondary-amine sites can differ from
older stored assignments without changing their standalone structures. No
bundled template or shipped example is renumbered by this policy change.

There are 55 definitions with warnings. Counts overlap: 47 have unspecified
tetrahedral stereo, eight have reviewed identity concerns, two lack some leaving
metadata, two retain chemistry declaration conflicts (`Hpic`, `Pca`), and one
has metadata referring to a missing site (`Hpic`). The two activation failures
also remain warnings. No configuration is invented for an unspecified center.

`library-curation.json` retains explicit review findings bound to an exact
definition hash. Changing a structure, its name, slots or chemistry metadata
detaches that finding until reviewed again. Regenerating the automatic audit
does not erase an applicable identity concern. Corrected records must be
reviewed against the stored finding before that record is removed or rebound.

## Runtime meaning

`pyPept.library_quality.quality_for_monomer(molecule)` performs a read-only
fingerprint lookup, without running activation or modifying any library. The
palette exposes that result as `quality`:

| Status | Meaning |
| --- | --- |
| `no_known_exception` | This exact definition matches the baseline and has no warning; informational notes can remain. |
| `review_required` | This exact definition matches the baseline and has an explicit warning to inspect. |
| `unreviewed` | This custom, renamed, changed, or otherwise unmatched definition has no applicable baseline. |

Only warnings appear as amber tile text. Hover previews retain informational
notes in neutral text and warnings separately. An unreviewed definition remains
usable. Its empty `issues` array must not be
read as evidence that its chemistry has been reviewed. A changed fingerprint cannot inherit the known
status of another record with the same label. Unsupported enhanced stereo also
cannot inherit a fingerprint produced by plain-SMILES normalization.

`audit_version` names the fingerprint/audit convention. Version 3 retains the full
name, severity, standalone preservation and curated findings, and records
explicit reactivation orientation choices under `canonical-sites-v1`.
The manifest records its
generating RDKit version for provenance. CI compares the per-record baseline on
the supported Python 3.9, 3.11 and 3.13 environments with their installed RDKit
versions. A future canonical-SMILES change that alters
a hash fails the audit regression and makes that runtime definition unreviewed
until checked. No cross-version scientific identity guarantee is implied.

## Reproducing and reviewing changes

Run the full audit into a candidate file before changing the packaged baseline:

```sh
python -m pyPept.library_quality --output /tmp/library-quality-candidate.json
python -m pytest tests/test_library_quality.py \
  tests/test_bundled_monomers.py::TestLibraryRoundTrip
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
TATA products, lysine sidechain acylation, and the corrected halide
substitution products. These controls do not prove the
identity of every imported record. Source attribution and scientific curation
of individual exceptions remain separate work.

## Recorded data concerns

The following entries retain explicit review warnings rather than guessed
structural corrections:

| Definition | Stored fact requiring review |
| --- | --- |
| `D_b3hTrp` | Its ring is hydrogenated despite the homotryptophan name. |
| `Hpic` | Its structure lacks the named carboxyl group. |
| `Pca` | Its carboxyl group is at position 3; the name specifies position 2. |
| `PhosOxScaffold` | Its carboxyl groups are meta to phosphorus; the name specifies ortho. |
| `TATA` | Unused arms restore as propionyl instead of the named acryloyl groups. |
| `bOH_Lys` | Its skeleton has five carbons rather than six. |
| `bOH_Val` | Its skeleton has six carbons rather than five. |
| `hArg` | Its full name and saturated diaminomethyl group differ from the homoarginine implied by the abbreviation. |

Name comparisons for pyridine carboxylates use [NCBI MeSH](https://www.ncbi.nlm.nih.gov/mesh/68010848)
and [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=98-98-6). Carbon-skeleton
comparisons use [PubChem 3-hydroxylysine](https://pubchem.ncbi.nlm.nih.gov/compound/3-Hydroxylysine)
and the [structural identification of beta-hydroxyvaline](https://pubmed.ncbi.nlm.nih.gov/7557542/).
These sources establish the named structures; they do not establish which
structure the original library author intended.

For `W`, the former declaration warning exposed a shared classifier bug:
aromatic nitrogen had been treated as a primary amine. Its R4 now resolves to
`aromatic_nh`. The builder rejects the unsupported R4-to-ac-R2 pairing and keeps
Connect disabled. The warning is absent because the classification is corrected.

The reported `ac-G-D.[R(4,1).G(2,1).D(2,1).am(2,1)]-S-am` retains its exact
assembled structure, 45 atoms and molecular weight 645.28. The former warning
counted peptide bonds inside brackets and proposed an unnecessary rewrite.
Brackets can contain a peptide branch; supported connections now use the shared
reaction policy. Browser checks confirmed the branch tab still highlights
exactly its four residues.

The audit covers the bundled definitions and their recorded issues. It does
not test every monomer pairing or establish synthesis feasibility.

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
