# CABILN notation

A CABILN string describes monomers and their connections. Each monomer is
resolved against the installed SDF library. A connection names attachment
slots, whose chemistry is checked against the reaction library.

## Backbone and caps

A dash always means **R2 on the left → R1 on the right**. Missing or occupied
slots are errors; a dash cannot mean a different bond.

```python
from pyPept.sequence import Sequence

Sequence('A-G-K')
Sequence('ac-A-G-K-am')
Sequence('ac-meA-G-am')  # preformed N-methylalanine
```

Caps have only the attachment sites defined in their library record.
N-methylalanine (`meA`) is an amino acid, not a methyl cap.

## Inline modification

`.Modifier(host_slot, modifier_slot)` attaches a modifier to the preceding
residue. Slot numbers refer to that monomer's library record.

```python
Sequence('fmoc-C.trt(4,2)-K.boc(4,2)-am')
```

Here Cys R4 attaches to trityl R2, and Lys R4 attaches to Boc R2.

## Crosslinks and cycles

The first endpoint names both slots. The second endpoint can use just the tag.
Each tag must occur exactly twice and must connect distinct monomers.

```python
Sequence('C.!1(4,4)-A-G-C.!1')
Sequence('ac-K.!1(4,4)-A-A-D.!1-am')
Sequence('!1-A-G-K-A-!1')
```

A terminal marker is a bond marker, not a residue. Adding residues to a cycle
must preserve which residues the two markers enclose.

## Sequential brackets

In `Host.[A(r,s).B(t,u)]`, A attaches to Host, then B attaches to A.
Each continuation starts from the preceding fragment.

```python
Sequence('K.[G(4,2).ac(1,2)]-A')
Sequence('ac-K.[E(4,4).AEEA(1,2).C20FA(1,2)]-G-am')
```

In the second example, standard glutamate `E` exposes its gamma carboxyl at R4.
`E_g` has a different assignment: its gamma carboxyl is R2 and its alpha
carboxyl is R4. These symbols cannot be substituted while keeping slot numbers.
The removed `gGlu` symbol should not be used.

An explicit child scope starts from the current monomer. Closing it restores
that monomer, so consecutive children are siblings. Children can contain their
own children. These two forms have identical connections:

```python
Sequence('K.[G(4,2).ac(1,2)]-A')
Sequence('K.[G(4,2).[ac(1,2)]]-A')
```

Here acetyl attaches to G in both forms. In contrast, these shapes differ:

```python
Sequence('K.[K(4,2).[K(1,2)].[K(4,2)]]-A')  # two children of the inner K
Sequence('K.[K(4,2).[K(1,2).[K(4,2)]]]-A')  # a child of a child
```

The older `[.arm]` and `[[...].X]` forms retain their existing meanings.
Changing continuation direction does not change which fragment is the source.
Reusing an occupied slot is invalid. See the
[canonical notation contract](canonical-notation-plan.md) for complete scope
rules and resource limits.

## Protected brackets and percent segments

`{...}` has the same assembly behavior as `[...]`, but layout formatting leaves
it in bracket form. The builder creates protected branches by default.

```python
Sequence('K.{G(4,2).ac(1,2)}-A')
```

`%` separates explicit chains. Every dash within each chain still means R2→R1.
Branches can attach through paired crosslink tags:

```python
Sequence('ac-K.!1(4,4)-G-am%C20FA-AEEA-E.!1')
```

This represents the same graph as the glutamate bracket example above.
The branch reads in reverse order because its backbone bonds have that direction.
Disconnected peptide chains can also be represented as separate segments.

`cabiln_to_branch()` converts a bracket only when a backbone chain can express
its connections. Single modifiers, protected brackets, and unsupported branch
shapes remain unchanged. `cabiln_to_bracket()` performs the reverse conversion.
These compatibility converters do not establish canonical graph ordering.

```python
from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch

bracket = 'E.[G(4,1).A(2,1)]-am'
branch = cabiln_to_branch(bracket)
restored = cabiln_to_bracket(branch)
```

## Output conventions

The requested output style is a layout preference. The modern
`pyPept.inputs.format_source(source, "bracket")` and `"percent"` paths check that
formatting preserves monomer occurrences, attachment slots, and the assembled
molecular structure. The library-free legacy text converters have separate
compatibility behavior.

For example, these forms describe the same peptide:

```text
K.[G(4,2).ac(1,2)]-A
K.!1(4,2)-A%ac-G.!1(2,4)
```

The bracket starts at G and continues through G's R1 to acetyl R2. The explicit
branch reads `ac-G`, since its dash must run from acetyl R2 to G R1. Each explicit
crosslink endpoint names its own slot first; its partner therefore reverses the
pair.

Formatting does not define a unique string for every molecular graph. Aliases
such as `Ala-Gly` and `A-G`, chosen crosslink tag names, and preserved source
groups can still yield different strings. Explicit chains have deterministic
ordering preferences; source groups also preserve editing intent. Protected
brackets remain protected. These are different guarantees from graph-wide
canonicalization or choosing one universal main chain.

Request canonical output explicitly to give equivalent resolved monomer graphs
the same spelling in each style:

```python
from pyPept.inputs import format_source

source = 'K.[G(4,2).ac(1,2)]-A'
assert format_source(source, 'bracket', canonical=True) == (
    'K.[G(4,2).[ac(1,2)]]-A'
)
assert format_source(source, 'percent', canonical=True) == (
    'K.!1(4,2)-A%ac-G.!1(2,4)'
)
assert format_source('Ala-Gly', 'bracket', canonical=True) == 'A-G'
```

Canonical bracket output writes every pendant scope explicitly. Canonical
percent output writes ordered R2→R1 chains and complete paired crosslink markers.
Both choose chain order, cycle cuts, sibling order, and labels from the resolved
graph. They preserve every monomer occurrence and numbered attachment.

Canonical export ignores source grouping and brace protection; it retains the
chosen monomer decomposition. Source-preserving layout formatting remains the
default. In the web app, select **Canonical** beside the notation buttons, then
choose `→[ ]` or `→%`. Undo restores the previous spelling.

The HTTP `/convert_notation` request accepts `canonical: true` and returns
convention, RDKit-version, and library/rule content identifiers with the result.
Canonical strings are qualified by that convention and binding. They are not
promised to remain byte-identical across unqualified RDKit or library upgrades.
See the [canonical contract](canonical-notation-plan.md) and its
[validation record](canonical-notation-validation.md).

The parser requires complete bracket entries and consistent crosslink slot pairs.
Library-free text converters preserve unsupported input unchanged. They do not
discard unknown annotations to make an input valid. See the
[normalization review](notation-normalization-review.md) for examples and checks.

## BILN and HELM

Declare legacy BILN explicitly to avoid ambiguous slot interpretation:

```python
Sequence('C(1,3)-A-A-C(1,3)', fmt='biln')
```

Use `pyPept.converter.Converter` for BILN/HELM interchange. The web app also
accepts these formats through its notation selector and conversion controls.

## Inline SMILES

The reverse converter can emit `<...>` fragments when a library token cannot
represent the source graph. Numbered dummy atoms in a full template name its
attachment slots. These fragments are checked by the same assembly pipeline.
A single fragment can contain more than one residue, so its chip cannot be
interpreted as a per-residue annotation.

## Validation

```python
report = Sequence.validate('C.!1(4,4)-A-G-C.!1')
assert report.ok, report.errors
```

Validation and successful assembly establish representability under the current
library. Library metadata, including specified stereochemistry, still determines
the resulting structure. Compare the exported molecule to a known reference
when the exact isomer or attachment position matters.
