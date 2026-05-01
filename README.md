# pyPept — CABILN Fork

**Chemistry Aware BILN (CABILN) extensions to pyPept**

> This is a development fork of [pyPept](https://github.com/Boehringer-Ingelheim/pyPept) by Boehringer Ingelheim.
> Original publication: [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2), *Journal of Cheminformatics*, 2023, 15:79.
> Original authors: Rodrigo Ochoa, J.B. Brown, Thomas Fox.
> Fork extensions: Cameron Beeley, 2026.

---

## What this fork adds

**CABILN** (Chemistry Aware BILN) is a superset of BILN that makes modified and cyclic peptide notation more expressive:

| Feature | Syntax | Example |
|---------|--------|---------|
| Inline cap / modification | `.Cap(host_r,cap_r)` | `fmoc-C.trt(4,1)-am` |
| Terminal cyclisation | `!n-...-!n` | `!1-A-A-A-A-!1` |
| Sequential bioconjugation bracket | `Res[.A(r,s).B(t,u)]` | `C[.Mal(4,1).DBCO(2,1)]` |
| Sidechain branch | `chain%%branch` | `K.!n(3,1)%%!n-PEG-am` |

In addition, a new SMARTS-driven monomer pre-activation pipeline replaces the manual CHUCKLES authoring workflow, and a CLI tool allows novel monomers to be registered directly from SMILES or CABILN.

---

## What still works from upstream pyPept

| Input | Status | Notes |
|-------|--------|-------|
| Linear BILN (`A-G-K`) | ✓ Works | Unchanged |
| Capped BILN (`fmoc-A-G-K-am`) | ✓ Works | Unchanged |
| FASTA (`PEPTIDE` → `P-E-P-T-I-D-E`) | ✓ Works | Via `run_pyPept --fasta` |
| HELM linear (`PEPTIDE1{...}$$$$V2.0`) | ✓ Works | Converter produces BILN |
| Old BILN crosslinks (`C(1,3)-A-C(1,3)`) | ✗ Rejected | Use `.!n(y,z)` CABILN notation |
| HELM crosslinks (`$PEPTIDE1,PEPTIDE1,...`) | ✗ Fails | Converter outputs old crosslink notation |

---

## Installation

```bash
pip install git+https://github.com/anagnorisis2peripeteia/pyPept.git
```

Or for development:

```bash
git clone https://github.com/anagnorisis2peripeteia/pyPept.git
cd pyPept
pip install -e ".[dev]"
```

Requires Python ≥ 3.9, RDKit, and BioPython.

---

## Quick start

### Linear and capped peptides (BILN)

```python
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit import Chem

seq = Sequence('fmoc-A-G-K-am')
mol = Molecule(seq)
print(Chem.MolToSmiles(mol.get_molecule(fmt='ROMol')))
```

### Inline sidechain modification (CABILN)

Attach a protecting group or cap to a sidechain R-group inline:

```python
# Cys with trityl protection on the thiol (R4)
seq = Sequence('fmoc-C.trt(4,1)-G-A-am')

# Lys with Boc on the ε-amine (R4)
seq = Sequence('fmoc-K.boc(4,1)-A-am')
```

### Head-to-tail cyclic peptides

```python
# Cyclic tetrapeptide — !1 marks the cyclisation endpoints
seq = Sequence('!1-A-G-K-A-!1')
```

### Sequential bioconjugation bracket

Attach multiple groups to a residue in sequence — each step's R-group refers to the preceding fragment:

```python
# Maleimide conjugation followed by DBCO on Cys
seq = Sequence('G-C[.Mal(4,1).DBCO(2,1)]-A')
```

### Sidechain branch

```python
# PEG branch on Lys ε-amine
seq = Sequence('A-K.!n(3,1)-G%%!n-PEG-am')
```

### FASTA input (CLI)

```bash
run_pyPept --fasta ACDEFGHIKLMNPQRSTVWY
```

### HELM input (linear only)

```bash
run_pyPept --helm 'PEPTIDE1{[ac].D.T.H.F.E.I.A.[am]}$$$$V2.0'
```

---

## Monomer pipeline

### Pre-activate a new monomer

Converts plain SMILES to a CHUCKLES fragment with auto-detected R-groups:

```python
from pyPept.interfaces.monomer_pipeline import pre_activate

result = pre_activate('N[C@@H](CS)C(=O)O')
print(result.chuckles)    # '[1*]N([3*])[C@@H](CS[4*])C([2*])=O'
print(result.leaving)     # {1: '[H]', 2: '[OH]', 3: '[H]', 4: '[H]'}
print(result.chem_types)  # {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod', 4: 'thiol'}
```

R1/R2 (backbone N and carboxyl C) are detected via a graph-topology rule that handles α-, β-, and γ-amino acids without hard-coded stereo assumptions. Sidechain slots (R3+) are assigned wherever chemistry permits.

### Register a new monomer via CLI

```bash
# From plain SMILES
pyPept-monomer-add --smiles "N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O" \
    --symbol Lys_Boc --name "Boc-Lysine"

# From CABILN (assemble first, then onboard)
pyPept-monomer-add --from-cabiln "K.boc(4,1)" --symbol Lys_Boc
```

The monomer is appended to `monomers.sdf` without overwriting existing entries.

### Programmatic registration

```python
from pyPept.interfaces.cli_monomer import register_monomer

result = register_monomer(
    smiles='N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O',
    symbol='Lys_Boc',
    name='Boc-Lysine',
)
```

---

## Crosslink notation (CABILN vs old BILN)

Old BILN used bare integer bond IDs: `C(1,3)-A-A-A-C(1,3)`. This is **rejected** in CABILN to avoid ambiguity with attachment point indices.

Use `.!n(host_r,cap_r)` instead:

```python
# Disulfide bridge between two Cys (R4 = thiol)
seq = Sequence('C.!1(4,4)-A-A-A-C.!1(4,4)')

# Lactam staple: Lys ε-amine (R4) to Asp sidechain carboxyl (R3)
seq = Sequence('K.!1(4,3)-A-A-A-D.!1(3,4)')
```

---

## CLI tools

| Command | Description |
|---------|-------------|
| `run_pyPept --biln <seq>` | Generate 2D/3D structure from BILN/CABILN |
| `run_pyPept --fasta <seq>` | Generate from FASTA single-letter codes |
| `run_pyPept --helm <seq>` | Generate from HELM (linear only) |
| `pyPept-BILN-validate --biln <seq>` | Validate a BILN/CABILN string |
| `pyPept-monomer-add --smiles <smi> --symbol <tok>` | Register monomer from SMILES |
| `pyPept-monomer-add --from-cabiln <seq> --symbol <tok>` | Register monomer from CABILN |

---

## Tests

```bash
pytest tests/test_bond_validation_and_assembly.py -v
```

260 tests covering bond validation, CABILN notation, monomer pipeline, and CLI.

---

## References

- [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2), *Journal of Cheminformatics*, 2023.
- [BILN — A Human-readable Line Notation for Complex Peptides](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00703), *J. Chem. Inf. Model.*, 2022.
