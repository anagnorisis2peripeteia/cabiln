# Repository notes

CABILN extends pyPept with numbered attachment sites, reaction-based assembly
and a web peptide builder. Python 3.9+ is supported. The browser uses JavaScript
without a UI framework.

## Commands

```bash
python -m pip install -e '.[dev,web]'
python -m pytest -m 'not distribution and not fuzz'
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
python tools/check_examples.py
cabiln
```

Run affected test modules first. The distribution check installs a fresh wheel
and needs package-index access. Browser setup is in
[tests/browser/README.md](tests/browser/README.md); generated tests and replay
are in [tools/README.md](tools/README.md#generated-tests).

## Code and chemistry

[CONTEXT.md](CONTEXT.md) defines the domain terms.
[Architecture](docs/architecture.md) maps module responsibilities.
[Notation](docs/notation.md) and [canonical output](docs/canonical-notation.md)
describe the format. [Test ownership](tests/AGENTS.md) identifies the suites.

- Monomer definitions supply registration, recognition, library tiles, attachment
  choices and assembly. Support new monomers through those definitions.
- Preserve declared R-group numbers. Amino-acid activation conventionally assigns
  R1 to backbone N, R2 to backbone C, R3 to the additional N attachment and R4+
  to sidechains. Resolve each definition's chemistry before connecting sites.
- A dash connects left R2 to right R1. Occurrence identity, atom index and site
  number are separate; repeated residues need separate identities.
- Preserve structure, charge, isotopes, specified stereo and atom ownership.
  Carboxyl and aldehyde sites can share an activated C(=O) graph; their leaving
  groups distinguish them. Exclude the backbone carboxyl from sidechain detection.
- Use the shared attachment, reaction and leaving-group implementations.
  Reaction routing and restoration rules live in `src/pyPept/data/`.
- Preserve the working builder's selection, panels, previews, highlights, history
  and exports. Test changed request lifetimes through their browser events.
- Keep library records and reaction data unchanged during structural refactors
  unless a demonstrated chemistry defect requires a separate correction.

## Maintenance

Use the repository's 88-column Python style, Black and isort's Black profile.
Keep Python 3.9 compatibility. Comments should explain chemistry or compatibility
constraints. Give shared rules one owner and keep public compatibility adapters
where callers still need them.

Update the relevant guide when behaviour, commands or ownership changes.
Keep dated measurements with their revision and environment.
Keep personal reviews, work plans and session reports outside the repository.
Do not publish those notes or link to copies in Git history.
Do not add AI co-author trailers to commits.
