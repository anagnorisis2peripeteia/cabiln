# Drawing and notation conversion performance

Interaction feedback and chemistry completion have separate budgets. Aim to
paint an acknowledgement within 200 ms; measure field responsiveness at the 75th
percentile separately for mobile and desktop, following the
[INP guidance](https://web.dev/articles/inp). A chemistry request may take longer
while leaving the controls responsive and the previous depiction visible.

The browser performance journey records both feedback and completed-result
timings for short edits/connections, a capped 60-residue chain and six nested
branches. Compare identical cases against a frozen revision with the same
browser, viewport and Python environment. Investigate a repeatable completion
regression above 20% before release; compare cold and warm work separately.
Do not infer production percentiles or chemistry throughput from first feedback.

Web rendering assembles the molecule without calculating a preliminary layout.
The default drawing uses Indigo through a detached, atom-mapped isomeric SMILES.
The returned structure must preserve the mapped graph and stereochemistry, and
provide a complete atom correspondence and valid two-dimensional coordinates.
Only those coordinates are copied to the original molecule; its atoms, bonds,
properties and monomer ownership remain authoritative. Unsupported structures
and failed validation use the existing CoordGen renderer.

Explicit odd layout seeds retain Indigo and even seeds retain the CoordGen
search. The first click of **Layout** selects seed 2, so it changes engines from
the default. The CoordGen search stops when it reaches zero atom overlaps: its
nonnegative score cannot improve further, and equal scores never replace the
current winner.

Structure verification assembles without coordinates because its comparison uses molecular structure.
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

When the source already has a current drawing, notation conversion retains its
SVG and MOL together. The formatter returns its validated mapping from target
occurrences to source occurrences, plus fresh target tabs and warnings. The UI
remaps residue ownership without changing atom indices. It retains the chosen
layout and zoom, clears old Build selections, and verifies the converted text
when Verify is active.

Reuse requires the same source, library and canonical context, a complete
occurrence mapping, an unchanged canvas, and the same accepted drawing at
response time. Edits, replacement drawings, external resizing and missing
metadata use the ordinary render path. Canvas comparison accounts for the space
taken by accepted tabs and status text. History still records notation changes;
Undo and Redo render their selected documents normally.

The formatter also assembles from the resolved source and target peptide objects
it already checked. This removes two repeated graph projections. Both chemical
assemblies, connection checks, monomer checks and exact product comparison remain.

The [notation reuse measurement record](../tools/benchmarks/results/notation-reuse-20261001.json)
compares the formatter with `60e9ff0`. Warm local medians across both notation
directions and both policies fell from 56–73 ms to 52–72 ms for Semaglutide,
and from 76–90 ms to 72–86 ms for Lixisenatide. Each median uses five requests.
All eight pre-existing response field sets match. The same record separates these local samples from the live
pre-change browser timings; it does not claim a hosted result from local data.

## Removal of duplicate work, 30 September 2026

These measurements compare `6c4ed23` with `9d2d946`. They are small local samples from 30 September 2026,
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

The comparison checked 42 complete HTTP responses and seven format/error
responses. HTTP regressions require one coordinate calculation per drawing and
none for verification. The source hashes, response checks and browser samples
are retained in the measurement record above.

## Faster default layout, 1 October 2026

The guarded Indigo default substantially reduces drawing time relative to
`9d2d946`. The same local HTTP capture clears the render cache before each call.
These are medians of three samples per revision on the machine described above.

| Operation | Previous default | Guarded Indigo default |
| --- | ---: | ---: |
| Semaglutide drawing | 729 ms | 56 ms |
| Lixisenatide drawing | 636 ms | 75 ms |
| Semaglutide reference, automatic format detection | 771 ms | 52 ms |

The machine had other active jobs; sample ranges and an earlier candidate capture are retained in the
[measurement record](../tools/benchmarks/results/fast-layout-20261001.json).
These local figures exclude network latency and hosted cold starts.

All 42 compared responses retain their nonvisual metadata and exported molecular
structure, including atom order, bond order and stereochemistry. Seven detection
and error cases also match. Drawing geometry changes. The engine-selection
experiment covered all fourteen examples and six additional chemical structures;
all twenty retained chemistry and ownership, with finite 2D coordinates and no
measured atom overlaps or bond crossings after normalization.

RDKit can add display-only hydrogens while preparing a drawing. These glyphs
inherit the original neighboring atom's SVG identity, so residue highlighting
includes them. Existing explicit hydrogens keep their own identity. Preparation
uses a separate drawing molecule; exported atoms and residue maps do not change.
The browser audit checks all 367 tabs across the fourteen examples, including
51 NH₂ groups, 759 stereo glyphs and both previously dim display hydrogens.

## Remaining costs

The alternate CoordGen layout can evaluate six additional candidates when its
initial drawing has overlaps. Semaglutide and Lixisenatide have nonzero raw
overlap scores, so the zero-score shortcut does not skip their initial search.
The local samples do not establish faster long-peptide CoordGen rerolls.
Every uncached drawing still needs its final coordinate calculation.
Superseding an active drawing still retires its process worker; an isolated 200-G cancellation
control needed about 3.25 seconds for replacement and admission retries.

Equivalent notation can change atom order, residue maps and the resulting SVG
when assembled afresh. Reuse therefore keeps the original SVG and MOL atom order
and uses the formatter's occurrence correspondence. Molecular equality alone
cannot supply that correspondence. Conversion of newly entered text still needs
its first target drawing, and all formatting retains two chemical assemblies.
