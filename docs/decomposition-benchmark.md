# Independent decomposition benchmark

Measured 2026-09-28. Fourteen independently specified structures, each tested as
supplied and after deterministic atom reordering, separate exact reconstruction
from useful residue recovery. This is a bounded diagnostic set, not a population
accuracy estimate or evidence that the search is globally optimal.

The frozen 28 September implementation passed all 28 cases with unchanged annotations.
The fixed baseline `f6c0423` recovered the expected partitions in 26; its
Sar/Gly decomposition exposed inconsistent backbone scoring, subsequently fixed
without a symbol-specific rule. Both preserve all specified chemistry. Datagrok's
comparator completed all calls faster, but changed connectivity on simple
independently reconstructable controls.

## Independent references

[Fixtures](../tools/benchmarks/decomposition_cases.json) contain literal source
SMILES, explicit atom ownership, attachment anchors and inter-residue bonds.
These annotations were authored from peptide chemistry and published annotations,
before running either importer. Neither CABILN assembly nor importer output
generated the reference molecules or partitions. The harness checks that the
reference partition covers every atom once and that its cuts are actual bonds.

The six published pairs were retrieved from PubChem's property and biologic
description APIs; their exact URLs, SMILES, HELM and InChIKeys are cached in
[pubchem_sources.json](../tools/benchmarks/pubchem_sources.json). No network access
is needed to run the benchmark.

| Reference | Independent expected partition |
| --- | --- |
| [Glutathione, CID 124886](https://pubchem.ncbi.nlm.nih.gov/compound/124886) | gamma-Glu, Cys, Gly; gamma-glutamyl bond |
| [Oxytocin, CID 439302](https://pubchem.ncbi.nlm.nih.gov/compound/439302) | Nine residues CYIQNCPLG, terminal amide cap, Cys1–Cys6 disulfide |
| [Oxidized glutathione, CID 65359](https://pubchem.ncbi.nlm.nih.gov/compound/65359) | Two gamma-Glu–Cys–Gly chains joined through cysteine sulfur |
| [Carnosine, CID 439224](https://pubchem.ncbi.nlm.nih.gov/compound/439224) | beta-Ala, His; preserve published imidazole tautomer |
| [N-acetylglycine, CID 10972](https://pubchem.ncbi.nlm.nih.gov/compound/10972) | Acetyl cap, Gly |
| [Leu-enkephalin, CID 461776](https://pubchem.ncbi.nlm.nih.gov/compound/461776) | Tyr, Gly, Gly, Phe, Leu, with specified stereo |

Eight further literal SMILES specify Ala–Gly, cyclic Ala–Gly, Sar–Gly, a
13C-O-methylserine between two glycines, alpha/epsilon-diglycyllysine,
unspecified-stereo Ala–Gly, two disconnected diglycines, and diglycine with an
additional registered whole-dipeptide monomer. The last case modifies only a
temporary library copy and accepts either two glycines or the registered
dipeptide. Molecular structure cannot uniquely establish a synthetic history.
Carnosine's exact His recognition is optional because tautomer-specific library
availability must not justify altering the published structure; its partition
and attachment anchors are still required.

## Scoring and limits

The [runner](../tools/benchmarks/decomposition.py) independently scores:

- Exact isomeric molecule equality, including isotope and charge. Inferred
  stereochemistry is separately disclosed and allowed only for the explicit
  unspecified-stereo fixture.
- Every original atom owned exactly once, and the expected residue atom sets
  recovered under a chirality-aware graph isomorphism.
- Agreement between reported source assignments and the actual emitted
  occurrences' atom ownership and library/synthetic status.
- Required R1/R2/sidechain source anchors, and all consumed attachment pairs
  matching the independently marked inter-residue bonds.
- Recognition flags for those expected residues. A whole-molecule synthetic
  escape cannot satisfy a multi-residue annotation.
- Conversion time, process wall time and peak RSS, separately.

Each case runs in a fresh process: 20 s wall, 15 s CPU soft/16 s hard, and 2 MiB
per output file. Timeouts terminate only that case's process group. Memory is
reported, not capped. Symmetry scoring is capped at 4096 isomorphisms and refuses
an exhausted comparison. One seed-17 atom permutation was measured; the runner
also supports seeds 41 and 73. No case hit a resource cap.

Eight scorer controls pass, including rejection of a correct whole-molecule escape,
correct atom partition with reversed peptide connectivity, and swapped attachment
anchors. Changed isotopes cannot hide behind allowed stereo inference; wrong
emitted occurrence indices also fail. An external success string without an
ownership witness stays unassessed for partition acceptance.

## Results from 28 September 2026

The archived source and library hashes remained unchanged throughout the run.
Python 3.11.15, RDKit 2026.03.6, macOS 27 arm64; each call includes cold importer
and library caches, with normal operating-system disk caching. Conversion time
excludes process startup and the harness's second reconstruction check. Baseline
SDF SHA-256 starts `bb0b9281dfa7`; complete hashes and per-case outputs are in
[baseline-f6c0423.json](../tools/benchmarks/results/baseline-f6c0423.json).
The final rework was copied from the frozen working tree into a separate source
root; its complete source hashes and unchanged-source check are in
[final-scoring-20260928.json](../tools/benchmarks/results/final-scoring-20260928.json).

| Measure, 28 calls | CABILN `f6c0423` | Final frozen rework | Datagrok `5b8151b` |
| --- | ---: | ---: | ---: |
| Returned a decomposition string | 28 | 28 | 28 |
| Independently exact molecule | 26 | 26 | 3 |
| Compatible only after unspecified stereo inference | 2, disclosed | 2, disclosed | 1, undisclosed |
| Incompatible rebuilt molecule | 0 | 0 | 14 |
| Reconstruction unavailable to verifier | 0 | 0 | 10 |
| Complete source-atom ownership | 28 | 28 | Unavailable |
| Ownership agrees with emitted occurrences | 28 | 28 | Unavailable |
| Expected residue partition and attachment semantics | 26 | 28 | Unavailable |
| Median / maximum conversion seconds | 0.424 / 3.037 | 0.308 / 0.611 | 0.096 / 0.171 |
| Median / maximum process wall seconds | 1.370 / 4.404 | 1.051 / 1.598 | 0.295 / 0.401 |

The [intermediate rework snapshot](../tools/benchmarks/results/rework-20260928.json)
is retained as history. Before the scoring correction it still recovered 26/28
expected partitions, with median conversion 0.290 s and maximum 0.558 s. This is
not the final implementation's result.

CABILN's largest measured RSS was 198 MB before and 210 MB after the rework.
Other local jobs ran during measurements, so observed runtime differences are
diagnostic rather than controlled speedup estimates. These timings also compare
different libraries and verification obligations; they are not a solver speed
ranking.
Datagrok used its bundled HELM library. It did not receive the CABILN-only
diglycine registration overlay. Complete external outputs are retained in
[mol-to-helm-5b8151b.json](../tools/benchmarks/results/mol-to-helm-5b8151b.json).

CABILN recovered every published partition in both atom orders, before and after
the rework. Its unknown isotope case returned a local synthetic residue with
R1/R2 and retained both known glycines. The two original annotation misses were
`CNCC(=O)NCC(=O)O`: both orders returned `Me_-G-G` rather than the expected
`meG-G`, despite the library containing `meG`. The molecule and atom coverage were
exact, but the residue granularity differed.

Review traced this to counting a methyl cap's R2–R1 connection as a backbone
bond even though the cap has no backbone. The corrected scoring and search bounds
require both endpoints to be backbone-capable; the existing larger-template
tie-break then selects `meG-G`. The final benchmark confirms the expected Sar/Gly
atom sets and R1/R2 anchors in both input orders. The reference annotations and
library data were not changed. Historical attachment assessment stays unavailable
where that reference partition was absent, rather than labeling alternate slots
invalid.

The minimal external counterexample is atom-order dependence:

```text
N[C@@H](C)C(=O)NCC(=O)O        -> PEPTIDE1{A.G}$$$$V2.0
C(O)(CNC([C@@H](N)C)=O)=O      -> PEPTIDE1{G.A}$$$$V2.0
```

The input molecules are identical. RDKit independently reconstructs the standard
HELM outputs and distinguishes Ala–Gly from Gly–Ala. Two disconnected diglycines
likewise become diglycine plus two free glycines. HELM containing custom monomers
outside RDKit's supported subset remains unassessed, not automatically incorrect.
The external public API supplies no ownership witness.

## Reproduction

First extract the fixed source into a new temporary directory, without switching
the working tree:

```bash
baseline=$(mktemp -d /tmp/cabiln-baseline.XXXXXX)
git archive df17590 | tar -x -C "$baseline"
.venv/bin/python tools/benchmarks/decomposition.py \
  --source-root "$baseline" --revision-label f6c0423 --permutations 1
.venv/bin/python tools/benchmarks/decomposition.py \
  --revision-label current-working-tree --permutations 1
.venv/bin/python tools/benchmarks/decomposition.py \
  --engine mol-to-helm --external-checkout /path/to/pinned/mol-to-helm \
  --permutations 1
.venv/bin/python -m pytest -q tests/test_decomposition_benchmark.py
```

Use the external checkout at
[`5b8151b4bab48891bd4403fe2e8ef644e65e3aed`](https://github.com/datagrok-ai/mol-to-helm/tree/5b8151b4bab48891bd4403fe2e8ef644e65e3aed)
with its requirements available; the adapter does not install it. The runner
prints its evidence directory and records failures in JSON; it exits zero after a
completed measurement, not only after perfect scores. `--case ID` selects a
discriminator and `--output DIRECTORY` pins the artifact location. Registration
overlays are local to those output directories; nonempty output directories are
refused. Archived converter responses were rescored after clarifying unavailable
attachment metrics and separately rebuilding their emitted ownership against the
pinned source and same library. Historical conversion timings were preserved;
each implementation snapshot has its own recorded measurements.

These results do not cover long therapeutic peptides, arbitrary non-peptide
chemistry, all reaction products, or library-order stability. Existing reaction
and large-peptide regressions remain separate obligations. They establish concrete
independent controls and a repeatable comparison surface, not a replacement
architecture decision.
