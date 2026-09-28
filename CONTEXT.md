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
