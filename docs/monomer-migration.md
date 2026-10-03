# Reprocessing an existing library

The bundled library uses `canonical-sites-v1` throughout. All 1,128 definitions
were reprocessed, including caps, scaffolds and noncanonical residues. Existing
site numbers changed on 124 definitions. Their standalone structures and the
reaction chemistry of every existing attachment were preserved.

Reprocessing restores each monomer, applies the current detector, combines its
findings with mandatory authored sites, and assigns canonical numbers. Stored
R1/R2 anchors retain the intended backbone or cap orientation. R3 is reserved
for a second attachment on the backbone nitrogen; sidechain sites use R4+.
Symmetry-equivalent handles are matched before combining capacities, so two
descriptions of the same maleimide handle do not create an extra site.

Handles outside automatic detection remain explicitly authored. Reprocessing
does not guess corrections to named chemical identities or unspecified stereo.
The [quality baseline](library-quality.md) records those remaining concerns.

## Migrate a custom library

Keep the original SDF and its companion `monomers.csv`. The command writes a
new directory and refuses to overwrite one that already exists:

```sh
python -m pyPept.monomer_migration \
  --source /path/to/original/monomers.sdf \
  --output /path/to/migrated
```

The output contains `monomers.sdf`, `monomers.csv` and `site-migration.json`.
The manifest binds the original and resulting files by hash and records every
old-to-new slot mapping. Names, types and CSV synonyms are retained. The CSV
also receives complete leaving groups, chemistry declarations and activation
policy. Reprocessing the resulting definitions again leaves them unchanged.

An explicit `build_library_from_csv(..., rebuild=True)` also reprocesses all
existing templates. Use the command above when migrating published notation:
its separate output and manifest retain the evidence needed to map old sites.
Ordinary library loading and CSV builds without `rebuild=True` do not renumber.

## Migrate CABILN references

Use the same original SDF and the generated output directory:

```sh
python -m pyPept.monomer_migration \
  --source /path/to/original/monomers.sdf \
  --output /path/to/migrated \
  --sequence 'meC.[C(3,4)]'
```

This prints `meC.[C(4,4)]`. The old R3 thiol is now R4. Other examples include
arginine's former R5 becoming R7, and TBAB's former R4/R5/R6 becoming R4/R7/R9.
Use the manifest for the complete mapping; do not apply a global number shift.

The migrator resolves the original notation against the original library,
maps connections by monomer occurrence, and writes from the resulting graph.
It assembles both versions and requires identical isomeric structures. Layout
can change when a connection is renumbered. A mismatch stops the migration.

For the previous bundled library, recover both original files from commit
`6fc8039d38bb0159be356a5ee612b8c2c8fd791e` in this repository. Custom libraries
need their own original files, including any registered monomers.

Saved projects retain their library context and reject incompatible definitions.
Recover each CABILN document or draft from the project JSON, migrate it with its
original library, then render, verify and save a new project. Bare notation has
no embedded library version: its numbers must be interpreted against the library
that produced it. Keep the original project until the migrated work is verified.
