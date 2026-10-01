# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

pyPept is a CABILN (Chemistry Aware BILN) fork of the Boehringer Ingelheim pyPept library. It converts single-string peptide notations into atomistic RDKit molecules. The fork adds: pinned R-group numbering (R1=backbone_n, R2=backbone_c, R3=backbone_n_mod, R4+=sidechain), SMIRKS-based bond assembly via reactions.yaml, inline cap/crosslink syntax, and a bundled monomer library with SMARTS-based auto-detection.

## Build & test commands

```bash
# Install (editable)
pip install -e ".[dev,web]"

# Run chemistry and HTTP regressions
python -m pytest -m "not distribution"

# Browser request-lifetime regressions
node --test tests/test_frontend.js

# Build and install a release artifact (requires package-index access)
python -m pytest tests/test_distribution.py

# Run a single test class or test
pytest tests/test_bond_validation_and_assembly.py::TestAssembly -v
pytest tests/test_notation_conversion.py::TestRoundTrips -v

# Verify every Sequence() call in README produces a valid molecule
python tools/check_examples.py

# Generate a candidate activation and metadata audit for every bundled monomer
python -m pyPept.library_quality --output /tmp/library-quality-candidate.json

# Start the live renderer web app (127.0.0.1:8732)
cabiln

# Register a new monomer from SMILES
pyPept-monomer-add --smiles "N[C@@H](CS)C(=O)O" --symbol Cys_check
```

## Architecture

See [docs/architecture.md](docs/architecture.md) for the web and conversion modules.
The former renderer now lives in `src/pyPept/web/`; reverse conversion lives in
`src/pyPept/smiles.py`. The `tools/live_renderer.py` file is a compatibility launcher.
Web registration is disabled by default; enable only on a trusted local instance.


### Core pipeline: string → molecule

```
CABILN string
  ↓  Sequence.__init__()            [sequence.py]
  │  ├─ _preprocess_cabiln()        segment splitting (%), terminal markers (!n)
  │  ├─ _expand_inline_caps()       .Cap(y,z) and .!n(y,z) → pendant chains + bond list
  │  ├─ _check_bond_chemistry()     validates each bond's element-pair chemistry
  │  └─ stores s_monomers[], s_bonds[]
  ↓
  ↓  Peptide.from_sequence()        [peptide.py]
  │  └─ resolved definitions, occurrence IDs and numbered endpoints
  ↓  Molecule(sequence_or_peptide)  [molecule.py]
  │  ├─ Sequence inputs adapt to the same resolved Peptide
  │  ├─ __assemble() assigns unique temporary labels to attachment dummies
  │  │   ├─ reaction_for_types() selects the shared reaction rule
  │  │   └─ run_bond_smirks() handles inter/intramolecular bonds
  │  └─ restore_leaving_groups() uses the same labels for unconsumed sites
  ↓
RDKit ROMol
```

### Key modules

| Module | Role |
|--------|------|
| `sequence.py` | CABILN/BILN parser, bond validation, monomer library loading from SDF |
| `molecule.py` | Assemble resolved peptide endpoints; retain Sequence compatibility |
| `interfaces/reaction_library.py` | YAML-driven reaction routing, `_CHEM_TYPE_REGISTRY` (SMARTS patterns), `infer_chem_type()`, `run_bond_smirks()` |
| `interfaces/monomer_pipeline.py` | `pre_activate()` (SMILES → CHUCKLES), `find_sidechain_slots()`, `build_library_from_csv()` |
| `interfaces/cli_monomer.py` | `register_monomer()` function and CLI entry point |
| `converter.py` | BILN ↔ HELM conversion (legacy, not part of CABILN pipeline) |

### Data files

| File | Format | Purpose |
|------|--------|---------|
| `data/monomers.sdf` | SDF with properties | Monomer library (CHUCKLES + leaving groups + chem_types) |
| `data/monomers.csv` | CSV | Authoring source for the core 52-monomer subset |
| `data/reactions.yaml` | YAML | 25 SMIRKS reactions (19 bond-forming + 6 terminal restoration) |
| `data/cap_reactions.yaml` | YAML | 100 cap-specific reactions (auto-applied) |

### CHUCKLES convention

Isotope-labelled dummy atoms encode attachment slots: `[1*]`=R1 (backbone N), `[2*]`=R2 (backbone C=O), `[3*]`=R3 (backbone-N mod), `[4*]`=R4+ (sidechain). The dummy's **neighbour** is the heavy atom where the bond forms. `_attachment_idx(mol, slot)` returns that neighbour; `_rgroup_atom_idx(mol, slot)` returns the dummy itself.

### SMIRKS reaction system

Reactions are defined in `reactions.yaml` with `reactant_pairs` that map `(chem_type_a, chem_type_b)` tuples to SMIRKS steps. The `REACTION_INDEX` dict is built at import time — adding a new reaction to the YAML file automatically makes it available without code changes.

Intramolecular ring closure uses RDKit's grouped-reactant syntax: `([A].[B]) >> [P]` called with `RunReactants((single_mol,))`. The isotope-swap fallback in `reaction_library.py` handles swapped dummy order. Targeted reaction labels explicitly match atomic-number-zero dummies, preserving isotopes on ordinary atoms.

### Monomer pre-activation

`pre_activate(smiles)` converts raw SMILES to CHUCKLES via:
1. **Backbone detection**: graph-topology shortest-path between amino N and carboxyl C
2. **R3 assignment**: if backbone N has ≥2 H (skipped for Pro)
3. **Sidechain detection**: `_CHEM_TYPE_REGISTRY` patterns in priority order, first-match-wins per atom
4. **Leaving group inference**: SMARTS-based rules determine what fragment (`[H]`, `[OH]`, `[Cl]`) to remove when placing the dummy

## Critical invariants

- **R-group numbering is pinned**: R1=backbone_n, R2=backbone_c, R3=backbone_n_mod, R4+=sidechain. Never varies by monomer. All CABILN notation, tests, and the web renderer depend on this.
- **Kekulize before dummy removal**: `molecule.py` Kekulizes the combined mol before removing any dummy atoms. This prevents aromatic-ring sanitization failures on His/Trp/etc. when an aromatic NH loses its dummy neighbour.
- **Backbone COOH excluded from sidechain scan**: After backbone detection, the backbone carboxyl's hydroxyl O is explicitly excluded so sidechain SMARTS won't re-match it (critical for Asp/Glu which have two COOH groups).
- **Carboxyl vs aldehyde disambiguation**: Both are `[CX3](=O)` in CHUCKLES. Distinguished by leaving group metadata: `[OH]` → carboxyl, `[H]` → aldehyde.
- **`take_largest: true`** in reactions.yaml: Filters out small byproduct fragments (NHS ring, N₂ from IEDDA). Required for any reaction with a leaving group.

## CABILN notation quick reference

```
-           backbone bond (amide by default)
.Cap(y,z)   inline cap attachment (host Ry to cap Rz)
.!n(y,z)    crosslink first endpoint (bond ID n)
.!n         crosslink second endpoint (inverse inferred)
%           segment separator (main chain first, branches after)
!n-...-!n   head-to-tail cyclisation (terminal markers)
[A(y,z).B(y,z)]  bracket multi-step conjugation
```

## Formatting

- Line length: 88 (flake8/black)
- Import sorting: isort with black profile
- No Co-Authored-By Claude tags in commits (causes CLA conflicts on OSS PRs)
