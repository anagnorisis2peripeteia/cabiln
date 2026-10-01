# Canonical CABILN

Canonical formatting gives equivalent resolved monomer graphs the same spelling
under a library binding and versioned convention. Request it with
`format_source(source, "bracket", canonical=True)` or the web app's **Canonical**
option. Default formatting preserves source layout.

The [notation guide](notation.md) covers ordinary input and conversion. This
reference defines recursive scopes, graph equality and the canonical writers.

## Recursive input rules

Each bracket scope has a host and a current monomer, initially its host.
Slot pairs mean `(host slot, new monomer slot)` for an attachment step.

1. The first monomer entry attaches to the scope's host and becomes its current monomer.
2. A flat `.Monomer(h,m)` continuation attaches to the current monomer and advances it.
3. A `.[...]` child attaches to the current monomer and opens a separate scope.
4. Closing that child restores the parent's current monomer.
5. Further child scopes at that position are siblings. Children can have children.
6. A crosslink marker belongs to the current monomer and does not advance it.

A marker-only scope therefore annotates its host. A marker's pair means
`(this endpoint's slot, partner endpoint's slot)`, wherever the marker appears.

Thus omitted brackets are accepted shorthand inside a sequential branch:

```text
Host.[A(h,a).B(p,b).C(q,c)]
Host.[A(h,a).[B(p,b).[C(q,c)]]]
```

Here the letters inside parentheses are placeholders for positive slot numbers.
These forms denote Host→A→B→C. They do not make B and C siblings.

Sibling branches and nested branches have different connections:

```text
K.[K(4,2).[K(1,2)].[K(4,2)]]-A
K.[K(4,2).[K(1,2).[K(4,2)]]]-A
```

```text
Sibling children:          A deeper child:

K0 ── A                   K0 ── A
└── K1                    └── K1
    ├── K2                    └── K2
    └── K3                        └── K3
```

These recursive inputs are accepted by `Sequence`. Missing closing delimiters
are errors; the reader must not guess
whether a missing delimiter changes the parent of a branch.

These rules also apply:

- Outside sequential brackets, `.Modifier(h,m)` suffixes independently modify
  the preceding backbone monomer. They do not form an implicit sequential arm.
- A dash always connects the left monomer's R2 to the right monomer's R1.
- The legacy `[.arm]` and `[[...].X]` spellings retain their existing meanings.
  Parse their compatibility productions into the same recursive records.
- Protected `{...}` groups retain their meaning and protection in layout mode.
- Treat `<SMILES>` as an opaque token while scanning notation delimiters.
  Square brackets and dots inside that token belong to SMILES.
- Consume the complete input. Validate both explicit inverse slot pairs,
  endpoint counts, slot existence, occupancy, and chemistry at every depth.

The parser uses explicit stacks, source spans and recorded parent/host ownership
to process nested scopes without relying on Python's call-stack limit.
The core limits input to 1,000,000 characters, nesting to 2,048 scopes, and one
scope parse to 100,000 entries. HTTP retains its 20,000-character request limit.

## What canonical equality means

Canonicalization operates on the **resolved monomer graph**. Each vertex is
one occurrence of a specific monomer definition. Each edge connects two numbered
attachment sites.
The library/reaction binding and canonical convention are part of the contract.

Equivalent graphs have a bijection preserving definitions, occurrences,
attachment slots, and connections. Source order, aliases, crosslink names,
source atom indices, and bracket grouping do not affect that equality.

Canonicalization preserves the chosen decomposition. A preformed modified
monomer and a base monomer plus modifier can assemble to the same molecule.
They remain different editable graphs. Formatting does not repeat recognition.

| Item | Canonical policy |
| --- | --- |
| Library aliases | Emit the canonical symbol of the resolved definition. |
| Distinct library records | Keep distinct unless explicitly declared aliases of one definition. |
| Synthetic monomers | Encode the full resolved template and attachment metadata canonically. |
| Repeated monomers | Keep every occurrence, including symmetric duplicates. |
| Opaque fragments | Keep each chosen fragment as one occurrence. |
| Slot labels | Preserve exactly; never renumber attachment sites for convenience. |
| Source tags | Replace with generated decimal tags. |
| Protected groups | Preserve in layout mode; ignore protection for plain-graph canonical export. |
| Chosen main chain | Preserve where required by layout mode; select from the graph for canonical export. |

Definition identity includes the resolved template, attachment chemistry,
terminal groups and role information used by the writers, including supported
stereochemistry, isotopes, charge and numbered dummies. Temporary synthetic
names and raw SMILES spelling cannot identify canonical definitions. Synthetic
templates remain distinct from named definitions unless explicitly resolved as
that definition; formatting does not perform a new library search.

`MonomerOccurrence.definition` carries a `ResolvedDefinition` with that identity.
Template snapshots are detached; canonical keys and tokens are computed lazily
and reused.

## Two canonical outputs

For the initial example, canonical bracket output is:

```text
K.[G(4,2).[ac(1,2)]]-A
```

Every pendant step receives its own explicit scope. Canonical output never
relies on omitted branch brackets. Input can still use the shorter spelling.

The corresponding canonical percent output is:

```text
K.!1(4,2)-A%ac-G.!1(2,4)
```

Percent output uses dashed R2→R1 chains and paired markers for other connections.
The chain reads `ac-G` because its dash connects acetyl R2 to glycine R1.
Both marker endpoints state their complete inverse slot pairs.

Both styles use canonical definition tokens, no optional whitespace, and decimal
slot/tag numbers without leading zeroes. Neither uses terminal-marker shorthand
or omitted crosslink slot pairs in canonical output. All supported monomer
tokens, including synthetic templates, must work inside recursive branches.

Bracket output still uses `%` between disconnected components and markers for
edges outside its spanning tree. Brackets alone cannot describe every cyclic
or multiply connected graph without revisiting occurrences.

## Ordering, paths, and cycles

Both writers use one graph-labeling implementation to produce an invariant
graph code and canonical vertex positions. Graph equality includes the chosen
definitions and numbered connections.

The graph code contains ordered definition identities and sorted edge endpoint
pairs `(vertex position, slot)`. The labeling process must resolve ties using
the whole graph. Input IDs, token spelling, recognition atom indices, and edge
insertion order cannot break output ties. Sorting symbols or local neighborhoods
alone is insufficient.

The writers apply these policies after labeling:

1. Canonicalize each connected component. Sort components by their graph codes;
   preserve the multiplicity of identical components.
2. Find maximal directed R2→R1 paths in each component, including isolated nodes.
   Slot occupancy ensures at most one predecessor and one successor per node.
3. For a pure directed cycle, start at its lowest canonical vertex position.
   Cut the incoming closing edge there and retain it as a residual connection.
4. Order paths by descending backbone-monomer count, then descending occurrence
   count, then their canonical vertex-position sequences.
5. Percent output writes every path in that order within its component.
   Every connection outside those dashed paths receives a paired tag.
6. Bracket output selects the first path as its primary dashed path and reserves
   all its vertices. Visit its vertices in path order to attach remaining nodes.
7. At each visited vertex, inspect non-primary edges in ascending order of
   `(own slot, neighbor canonical position, neighbor slot)`.
   First discovery creates an explicit child scope; recursively visit that child.
   Already discovered vertices receive no duplicate occurrence.
8. Every edge outside the selected primary path and child tree receives a tag.
   This includes ring closures, links between branches, and extra parallel edges.
9. Sort residual edges by component position and their ordered canonical
   endpoint pairs. Number them `!1`, `!2`, and so on across the document.
10. Emit a monomer's markers in tag-number order, before its child scopes.
    Emit child scopes in discovery order, then resume the parent path.

The tree is selected before writing text or assigning residual tags. The two
styles can use different tags for the same connection because their residual
edge sets differ. Tag numbers are serialization labels, not persistent bond IDs.

Symmetric occurrences can share a unique canonical string without possessing a
unique source-to-canonical correspondence. Each emission must return a valid
bijection to its original occurrences. Identical disconnected components obey
the same rule: exchanged source correspondences cannot change the output text.

### Labeling implementation and versioning

The implementation uses RDKit for a colored auxiliary
graph. Monomer vertices carry full definition identities. Each connection uses
two port vertices carrying the respective slot labels. Distinct label prefixes
separate monomer and port vertices. This retains parallel links and endpoint
orientation without treating residue connections as atomic chemistry.

RDKit's fragment ranking accepts custom atom symbols. Its
[canonical ordering example](https://www.rdkit.org/docs/Cookbook.html#reorder-atoms)
and [ranking interface](https://rdkit.org/docs/cppapi/namespaceRDKit_1_1Canon.html)
describe those interfaces. `canonical.py` adapts them to monomer graphs.

`canonical_convention()` identifies `cabiln-graph-v1`,
`rdkit-colored-port-graph-v1`, and the complete RDKit version. Canonical HTTP
responses also include content hashes for the selected monomer library, aliases,
reaction rules, and cap rules. Hashes are cached by file revision and expose no
local paths.
RDKit documents
[changes to canonicalization across releases](https://rdkit.org/docs/BackwardsIncompatibleChanges.html).
Ordering can change after dependency upgrades. Check stored canonical fixtures
and version any change to the output convention.

## Integration

`notation.py` reads recursive syntax and `Sequence` resolves its monomers and
connections. `canonical.py` labels the resulting `Peptide` and writes the selected
style. `inputs.format_source` reparses and verifies the emitted definitions,
connections and assembled structure.

`Serialization(text, occurrence_order, layout)` supplies new source spans and
a mapping back to the original occurrences. Canonical ranks are not persistent
occurrence IDs. Symmetric occurrences can have more than one valid mapping.

The core calls percent style `percent`; the HTTP adapter calls it `branch`.
Library-free legacy converters retain their permissive compatibility behaviour;
canonical export either produces verified output or reports an error.

## Regression properties

The following properties apply to both styles:

- **Equivalence:** implicit/explicit branch forms, aliases, renamed tags,
  permuted siblings, and reordered segments converge to the same string.
- **Distinction:** sibling/deeper graphs and changed attachment slots stay
  different. Different selected decompositions also stay different.
- **Idempotence:** canonicalizing canonical output does not change its bytes.
- **Round trips:** both writers preserve every occurrence, definition, and
  endpoint pair; conversion through either style reaches the same canonical
  result in the other. Assembled structures also compare exactly.
- **General graphs:** cycles, extra links, repeated monomers, and disconnected
  duplicates serialize once per occurrence with no lost connections.
- **Synthetic chemistry:** equivalent full templates converge while specified
  stereo, isotopes, charge, numbered sites, and terminal groups remain intact.
- **Editing:** selecting either of two identical monomers or caps edits only
  that occurrence. Closing nested branches restores the correct parent.
- **Compatibility:** accepted legacy syntax, BILN/HELM adapters, inline modifiers,
  protected groups, and molecular brackets inside `<SMILES>` keep their meaning.

Literal expected outputs, graph permutations and independently assembled
structures are covered in [test_canonical_notation.py](../tests/test_canonical_notation.py)
and the notation, editor and browser suites.
