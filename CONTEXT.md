# Peptide construction and recognition

CABILN describes peptides and related structures through monomers and their
connections. A molecular structure can admit several valid descriptions.

## Language

**Monomer definition**:
A named building block with a defined molecular structure, attachment sites,
and terminal groups. Several occurrences can use the same definition.

**Monomer occurrence**:
One particular building block in a peptide. Two glycines in the same peptide
are separate occurrences even though they share a definition.

**Attachment site**:
A numbered connection point on a monomer definition. Its chemical role and
terminal group determine how it participates in supported connections.

**Connection**:
A link between two attachment sites on monomer occurrences. Supported reaction
chemistry determines the resulting atoms and bonds.

**Atom ownership**:
The association of a product atom with a monomer occurrence. Ownership supports
highlighting and editing; it does not prove a unique historical synthesis.

**Decomposition**:
An interpretation of a molecular structure as monomer occurrences and their
connections, including explicitly identified unknown regions.

**Unknown region**:
A connected part of the source structure without a supported library match.
An editable unknown region has justified attachment sites; an opaque encoding
preserves structure without claiming those monomer boundaries.

**Notation layout**:
The choice of chains, branches, bracket grouping and labels used to write a
peptide description. Equivalent layouts preserve the occurrences and connections.

**Branch host**:
The monomer occurrence from which an attached branch starts. Two sibling
branches have the same host; a child branch starts from an occurrence inside
another branch.

**Bracket scope**:
The part of a notation describing an attached branch and its descendants.
Closing the scope returns to its host; a scope can contain several child scopes.

**Implicit continuation**:
A branch step whose own brackets are omitted in the written notation. Its host
and attachment sites remain the same as in the equivalent explicit nesting.

**Resolved monomer graph**:
A peptide description with specific monomer definitions, distinct occurrences,
and numbered-site connections. Different decompositions of the same molecule
can give different resolved monomer graphs.

**Canonical CABILN**:
The unique spelling of a resolved monomer graph within one output notation,
library binding, and versioned convention. It does not identify one unique
decomposition of every molecular structure.

**Library quality baseline**:
Versioned compatibility observations tied to an exact definition fingerprint.
Known exceptions remain explicit; a new or changed definition is unreviewed.
The baseline does not certify chemical identity or biological activity.

**Project context**:
The library, reaction rules and canonical convention under which source and
recognition evidence were resolved. Each document, notation draft and reference
retains its own context. Saved resolution signatures permit unrelated library
additions without silently accepting changed selected definitions.
