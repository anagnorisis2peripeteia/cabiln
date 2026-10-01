# Molecular structures and monomer recognition

The converter finds one verified interpretation of a molecular structure using
the selected monomer library. Newly registered monomers enter recognition,
palette discovery, assembly, and editing through that same library.

A lossless molecular encoding and a useful monomer decomposition are separate
requirements. Converting an entire peptide into one synthetic token preserves
its structure but does not recover editable residues. Tests check both.

## Algorithm

1. Parse every input component without changing its charge, isotopes, tautomer,
   or specified stereochemistry. Canonical atom order makes selection independent
   of SMILES traversal; an explicit map retains original input atom indices.
2. Filter templates by a necessary heavy-atom subgraph, then compile viable
   templates in their possible attachment states. The filter checks only
   elements and connectivity; it cannot reject a state for charge, isotope,
   stereo, hydrogen count, or bond order. Consumed slots match external
   neighbors; unused slots receive the same leaving groups as assembly. Exact
   matching checks the full chemistry and retains owned atoms, numbered slots,
   anchors, and effective attachment chemistry. No list of preferred monomer
   names is needed to ingest a new monomer.
3. Search for disjoint, compatible candidate occurrences across the entire
   component. Known neighbors must agree on their shared boundary. Search
   prioritizes known atom coverage and R2-to-R1 connections between occurrences
   with supported backbone chemistry, with stable ordering for equivalent
   alternatives. A methyl cap's R2 alkylation is not a peptide backbone bond.
   Counting it that way split sarcosine into methyl plus glycine despite an exact
   registered sarcosine template; the independent benchmark exposed this error.
4. Preserve unmatched connected regions with source attachment boundaries.
   Free termini require an unambiguous amine/carboxyl pair under the same
   topology convention used by ingestion. Unsupported boundaries reject the
   interpretation. They do not justify an invented residue or attachment.
5. Construct the shared occurrence/site/connection model. Its serializer lays out
   chains, branches, cycles and scaffolds, returning the output occurrence order.
   Independently assemble the notation and compare it with the input. Only an
   accepted interpretation can set the search's best coverage. This verification
   reuses versioned library definitions and does not generate drawing coordinates.
6. For supported reactions that change a monomer core, consider explicit
   precursor proposals with atom provenance, then use the same search and
   forward verification. Current adapters cover CuAAC, SPAAC, and ring-closing
   metathesis. Virtual precursor/byproduct atoms never become claimed input
   atoms.

Checking admissibility during search matters. For example, an isolated nitrogen
inside an unknown histidine tautomer can match a tiny amide template. Maximizing
the number of matched atoms before checking the remaining region would favor
that misleading cut and discard the useful whole-residue interpretation.

Cheap proposal search and admitted ownership search retain distinct traversal
and budget semantics. Both use the same reciprocal-port, unknown-atom and
cover-retention rules. Candidates and verified notation decisions can be reused
for one active proposal within a conversion; switching proposals releases that
state. There is no cache across input molecules or library revisions. The
separate budgets preserve partial progress at tight search limits.

Some local failures also provide a safe search constraint. Current synthetic
residues require a carbonyl attachment site. A sealed unmatched region without
one cannot be emitted this way. The search checks whether other library
candidates can explain that region before reusing the failure to reject its
surrounding ownership choices. A generic assembly error does not justify this
pruning. This distinction matters for large structures with many equivalent
choices far from one unknown residue.

## Result contract

`convert_smiles` returns a `ConversionResult`:

- `recognition_status="complete"`: known library templates cover every source
  atom in the selected interpretation.
- `recognition_status="partial"`: one or more local synthetic regions preserve
  source chemistry and supported boundaries alongside any known occurrences.
- `recognition_status="unresolved"`: at least one component needs an opaque
  encoding. Such a component has no claimed residue assignments.
- `assignments`: each recognized or local synthetic occurrence, its original
  input atom indices, its source attachment anchors, and its explicit
  `residue_index` in the returned CABILN. Indices remain correct when another
  component is unresolved or HTTP output uses bracket notation.
- `details`: the legacy primary-backbone view. It is not a list of all branches,
  caps, or scaffold occurrences; consumers needing all occurrences use
  `assignments`.
- `search_complete`: whether configured search/proposal limits were reached.
  A verified result can have an incomplete search. Diagnostics disclose this.
- `inferred_stereo`: whether the selected library supplies stereochemistry
  omitted by the input. Every specified input configuration must survive.

These fields do not establish a unique synthesis history. For example, an
ethylamide terminus can also match a reduced glycine template with a free
attachment slot. That is an alternative library interpretation with different
future editing capabilities, not merely a different spelling. Symmetric
scaffold arms and equivalent aliases also permit different valid notation.

The renderer derives tab ownership from the assembled occurrences. Reaction
atom maps preserve ownership through transformed junctions. Conversion source
assignments and renderer maps are checked against one another, including click
reactions and notation that reorders branch occurrences.
`convert_smiles(..., notation="bracket")` selects the layout inside serialization,
so assignment indices do not require a later molecular-isomorphism search.

## Limits

Recognition uses explicit pattern, match, candidate, search-state, solution,
slot-assignment, and retained-alternative budgets. Reverse proposals have their
own limits. Atom ownership is searched separately from template/slot choices
that share the same ownership. This prevents equivalent aliases from delaying
the exploration of a different, useful residue partition.
Limits prevent unbounded alias and symmetry searches; they do not prove that
all interpretations were examined. Unsupported transformation chemistry can
remain unrecognized even when forward assembly supports it. An intramonomer
imide closure is a recorded example: the tested structure retains its known
caps around a local synthetic residue.

Compiled-state caches retain at most two subsets and 50,000 query atoms in
total. A larger current input can require more temporary memory. Cache keys
include the selected library version, so registration invalidates stale data.

Opaque fallback is accepted only when full forward assembly preserves the
molecule. Its synthetic slots are encoding devices, not discovered residue
boundaries. If no exact representation is available, conversion fails with a
clear error. Unsupported counterions and enhanced AND/OR stereo groups are
never silently discarded.

The selected library is checked before and after conversion. A concurrent
library update invalidates the attempt and asks the caller to retry. Existing
library metadata issues listed in [library quality](library-quality.md) remain
separate data work.

## Verification and measurements

The [independent benchmark](decomposition-benchmark.md) checks external reference
structures, expected monomer ownership, attachment sites and atom-order changes.
The chemistry and browser suites also cover registration, unknown regions,
budget exhaustion and editable large-peptide imports.

The necessary-core filter was checked against unfiltered compilation for
mixed backbones, cycles, branches, unknown residues, CuAAC and metathesis
proposals. Every exact candidate was retained in those cases. The fresh-process
comparison recorded at `f6c0423` on 28 September 2026 used
`NCC(=O)N[C@@H](C)C(=O)O`. Filtering reduced peak resident memory from
636.3 MiB to 144.7 MiB, and call time from 4.76 to 1.27 seconds. The converter
compiled 128 viable states instead of 14,618 whole-library states. These are
measurements on the review machine, not memory or latency guarantees.
