# Drawing and notation conversion performance

Web rendering assembles the molecule without calculating a preliminary layout.
The drawing routine supplies the final coordinates. Structure verification
assembles without coordinates because its comparison uses molecular structure.
Reference format detection keeps its existing precedence and validates assembly;
the web caller opts out of its preliminary layout. Other callers retain the
existing depiction default.

When a user converts notation before the typing timer fires, the formatter
validates the entered source and the UI draws the converted result once. A failed
conversion still draws the entered source and retains the conversion error.
Active calculations finish before conversion, and newer actions invalidate stale
results. History records the exact text entered before an immediate conversion.
Undo's ordinary renderer can subsequently normalize legacy spelling. Normalization
by an already active drawing keeps its existing history behavior.

## Measurements

The baseline is `6c4ed23`. These are small local samples from 30 September 2026,
using an Apple M4 Pro, Python 3.11.15 and RDKit 2026.03.6. Independent baseline and
candidate processes use the same inputs and environment. Render misses clear the
existing cache; identical request hits are measured separately.

| Operation | Baseline | Candidate |
| --- | ---: | ---: |
| Semaglutide drawing | 1,104 ms | 581 ms |
| Lixisenatide drawing | 855 ms | 446 ms |
| Cyclic lipid drawing | 29 ms | 18 ms |
| Semaglutide reference, automatic format detection | 1,150 ms | 580 ms |
| Semaglutide exact verification | 560 ms | 29 ms |

These are medians of three calls per revision. Short G/A-G requests varied by a
few milliseconds and do not establish an improvement. The default Semaglutide
profile spent about 525 ms on preliminary coordinates and 532 ms on final
coordinates, while chemical assembly took about 14 ms.

The initial browser profile covered 54 operations in an isolated process-mode
server. Response-to-DOM work had a median of 0.8 ms and maximum of 2.9 ms, with
no tasks exceeding 50 ms. First notation conversions of long peptides spent most
of their time in the subsequent redraw. Clicking conversion during the typing
delay added an intermediate redraw as well.

| Browser event to completed drawing and two frames | Baseline | Candidate |
| --- | ---: | ---: |
| Changed Semaglutide input | 1,353 ms | 803 ms |
| First Semaglutide conversion to percent notation | 1,197 ms | 681 ms |
| Semaglutide conversion immediately after typing | 2,596 ms | 668 ms |
| Changed Lixisenatide input | 1,086 ms | 664 ms |

These are individual operations from complete 54-operation browser runs, with
all displayed structures checked independently after timing. Both layout and
canonical notation policies were exercised in both directions. Two canonical
Semaglutide operations in the candidate run were slower outliers; the complete
record retains them alongside the faster operations. Five subsequent repetitions
of each exact request did not reproduce those outliers: baseline rendering took
1.12–1.26 seconds and candidate rendering took 0.58–0.65 seconds, with identical
responses. The cause of the two original outliers was not isolated.

A frontend-only control on the old backend reduced immediate conversion from
2,596 ms to 1,270 ms. It removed the intermediate drawing; removing duplicate
coordinates then reduced the combined result to 668 ms.

[Recorded samples, source hashes and response checks](../tools/benchmarks/results/drawing-conversion-20260930.json)
distinguish backend calculations from browser event-to-paint measurements.
Local timings exclude hosted cold starts and network latency.

## Preservation checks

All 42 compared request responses are identical across revisions. The comparison
includes every built-in example at the default layout, five representative
peptides at both alternate layout seeds, explicit and automatic reference
formats, and exact, unspecified-stereo and mismatching verification cases.
It compares the complete response, including SVG, MOL export, residue ownership,
branch tabs, warnings and context. Seven format-detection and error responses
also remain identical.

The HTTP regression counts actual coordinate calculations: one for each main or
reference drawing and none for verification. All four cases fail on the baseline
and pass with this change. Frontend contracts cover request order, active
normalization, Undo, failure fallback, replacement clicks, reference updates,
Verify and edits that preserve the trimmed source. The latter prevents a
cancellation callback from scheduling two competing drawings.

Notation conversion retains occurrence, connection, monomer-definition and
assembled-structure checks. Library bindings, recognition and ingestion are
unchanged. No additional cache or runtime dependency is introduced.

## Remaining costs

The alternate CoordGen layout still evaluates six candidates. The default long
peptide layout still needs its final coordinate calculation. Superseding an
active drawing still retires its process worker; an isolated 200-G cancellation
control needed about 3.25 seconds for replacement and admission retries.

Equivalent notation can change atom order, residue maps and the resulting SVG.
Semaglutide's percent form did so in the comparison. Reusing a drawing based only
on canonical molecular SMILES would require an explicit correspondence between
atom indices and monomer occurrences. Molecular equality alone does not provide
that correspondence.
