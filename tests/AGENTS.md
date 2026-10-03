# Test ownership

Keep one primary test for each contract. Another layer needs a distinct failure
mode, such as HTTP field mapping, persistence, cancellation, or wheel packaging.

The chemistry suites are organized by responsibility:

| File | Responsibility |
| --- | --- |
| `test_bond_validation_and_assembly.py` | Bond validation, assembly products, and errors |
| `test_monomer_activation.py` | Raw input, activation, and CSV ingestion |
| `test_monomer_migration.py` | Authored-site reprocessing, numbering migration, and preserved products |
| `test_attachment_reactions.py` | Attachment inference and reaction-family products |
| `test_leaving_group_restoration.py` | Restoring leaving groups on activated monomers |
| `test_sequence_parsing.py` | Legacy CABILN lowering and Sequence validation |
| `test_notation_conversion.py` | Alternate notation forms and round trips |
| `test_monomer_cli.py` | Monomer CLI registration and argument handling |
| `test_smiles_to_cabiln.py` | Reverse recognition and editable atom partitions |
| `test_bundled_monomers.py` | Bundled chemistry, aliases, and library audit |

`_chemistry_oracles.py` contains the reused assembly and atom-partition checks.
Keep fixture tables and specialized assertions in their owning test modules.
Run these suites with pytest; there is no separate script runner.

- Raw monomer activation belongs to `test_monomer_activation.py`'s
  `TestMonomerPreActivate`; restoration belongs to
  `test_leaving_group_restoration.py`'s `TestRestoreRgroups`. Add labeled literal
  rows for new molecules before adding another copy of setup and assertions.
- Retain independent chemistry expectations: attachment atoms, leaving groups,
  products, stereo, charge, isotope, and atom ownership. The library audit does
  not replace these checks. Read helpers before calling a test assertion-free.
- Compare exact input and expected-fact sets before consolidating cases. Keep
  distinct raw spellings, D/L forms, branch directions, and error paths. Carry
  extra assertions before deleting a weaker test. Use ordinary pytest tables;
  avoid generated test functions or a custom scenario interpreter.
- A converter test that checks both literal forms and both directions already
  owns that pair. Preserve molecular equivalence between different forms;
  assembling two strings already asserted equal supplies no independent oracle.
- Exercise shipped rules, not copies of their SMARTS or reaction definitions.
  Keep independently authored positive and negative molecules as expectations.
- GUI races belong to their actual event/request lifecycle. Core chemistry
  tests do not replace those tests or the HTTP, storage, and worker boundaries.
- `test_runtime.py`, `test_launch_chemistry.py`, and `test_projects.py` also run
  against the installed release wheel. Preserve that CI selection when moving
  tests, including direct AND/OR stereo rejection in the launch suite.

Run the affected owners, then the normal Python and Node gates. Use independent
preservation review and focused fault controls for assertions moved during a
large consolidation. Report removed cases separately from retained table rows.
