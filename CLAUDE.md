# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

pyPept is a CABILN (Chemistry Aware BILN) fork of the Boehringer Ingelheim pyPept library. It converts single-string peptide notations into atomistic RDKit molecules. The fork adds: pinned R-group numbering (R1=backbone_n, R2=backbone_c, R3=backbone_n_mod, R4+=sidechain), SMIRKS-based bond assembly via reactions.yaml, inline cap/crosslink syntax, and a 1001-monomer library with SMARTS-based auto-detection.

## Build & test commands

```bash
# Install (editable)
pip install -e ".[dev]"

# Run the main test suite (572 tests, ~70s)
pytest tests/test_bond_validation_and_assembly.py -v

# Run a single test class or test
pytest tests/test_bond_validation_and_assembly.py::TestDisulfide -v
pytest tests/test_bond_validation_and_assembly.py::TestRoundTrips::test_round_trip_library -v

# Verify every Sequence() call in README produces a valid molecule
python tools/_test_readme_examples.py

# Full library round-trip (1001 monomers → CABILN → Molecule → SMILES)
python tools/full_library_roundtrip.py

# Start the live renderer web app (0.0.0.0:8732)
python tools/live_renderer.py

# Register a new monomer from SMILES
pyPept-monomer-add --smiles "N[C@@H](CS)C(=O)O" --symbol Cys_check
```

## Architecture

### Core pipeline: string → molecule

```
CABILN string
  ↓  Sequence.__init__()            [sequence.py]
  │  ├─ _preprocess_cabiln()        segment splitting (%), terminal markers (!n)
  │  ├─ _expand_inline_caps()       .Cap(y,z) and .!n(y,z) → pendant chains + bond list
  │  ├─ _check_bond_chemistry()     validates each bond's element-pair chemistry
  │  └─ stores s_monomers[], s_bonds[]
  ↓
  ↓  Molecule.__init__(sequence)    [molecule.py]
  │  ├─ __combine_all_monomers()    CombineMols into one RWMol
  │  ├─ __add_bonds_to_mol()        SMIRKS reactions via reaction_library
  │  │   ├─ relabel all dummies to globally unique isotopes
  │  │   ├─ infer_chem_type() for each attachment atom
  │  │   ├─ REACTION_INDEX[(type_a, type_b)] → YAML entry
  │  │   └─ run_bond_smirks() with inter/intramolecular handling
  │  └─ __restore_and_remove_rgroups()  cap unbound sites, Kekulize-first
  ↓
RDKit ROMol
```

### Key modules

| Module | Role |
|--------|------|
| `sequence.py` | CABILN/BILN parser, bond validation, monomer library loading from SDF |
| `molecule.py` | RDKit molecule assembly via SMIRKS reactions, R-group restoration |
| `interfaces/reaction_library.py` | YAML-driven reaction routing, `_CHEM_TYPE_REGISTRY` (SMARTS patterns), `infer_chem_type()`, `run_bond_smirks()` |
| `interfaces/monomer_pipeline.py` | `pre_activate()` (SMILES → CHUCKLES), `find_sidechain_slots()`, `build_library_from_csv()` |
| `interfaces/cli_monomer.py` | `register_monomer()` function and CLI entry point |
| `converter.py` | BILN ↔ HELM conversion (legacy, not part of CABILN pipeline) |

### Data files

| File | Format | Purpose |
|------|--------|---------|
| `data/monomers.sdf` | SDF with properties | 1001-monomer library (CHUCKLES + leaving groups + chem_types) |
| `data/monomers.csv` | CSV | Authoring source for the core 52-monomer subset |
| `data/reactions.yaml` | YAML | 25 SMIRKS reactions (19 bond-forming + 6 terminal restoration) |
| `data/cap_reactions.yaml` | YAML | 100 cap-specific reactions (auto-applied) |

### CHUCKLES convention

Isotope-labelled dummy atoms encode attachment slots: `[1*]`=R1 (backbone N), `[2*]`=R2 (backbone C=O), `[3*]`=R3 (backbone-N mod), `[4*]`=R4+ (sidechain). The dummy's **neighbour** is the heavy atom where the bond forms. `_attachment_idx(mol, slot)` returns that neighbour; `_rgroup_atom_idx(mol, slot)` returns the dummy itself.

### SMIRKS reaction system

Reactions are defined in `reactions.yaml` with `reactant_pairs` that map `(chem_type_a, chem_type_b)` tuples to SMIRKS steps. The `REACTION_INDEX` dict is built at import time — adding a new reaction to the YAML file automatically makes it available without code changes.

Intramolecular ring closure uses RDKit's grouped-reactant syntax: `([A].[B]) >> [P]` called with `RunReactants((single_mol,))`. The isotope-swap fallback (lines ~273-285 of reaction_library.py) handles cases where the two dummies appear in swapped order within the assembled molecule.

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
