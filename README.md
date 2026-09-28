# CABILN / pyPept

CABILN is a peptide notation and interactive builder built on
[Boehringer Ingelheim's pyPept](https://github.com/Boehringer-Ingelheim/pyPept).
Its monomer library supplies the builder's vocabulary. Ingesting a structure
detects its attachment sites and chemistry, creates its library tile, and makes
those sites available for building. New monomers with supported attachment
chemistry need no code specific to their names.

CABILN describes backbones, caps, branches, and crosslinks. The assembler uses
the monomer structures and reaction definitions to produce an RDKit molecule.

The web app provides structure drawing, residue selection, bond editing,
SMILES conversion, structure comparison, and PNG/MOL export.
The hosted application is at [cabiln.onrender.com](https://cabiln.onrender.com/).
A local checkout can contain changes that have not been deployed there.

## Install and run

Python 3.9 or newer is required. From a checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[web]'
cabiln
```

Open `http://127.0.0.1:8732`. On Windows, activate with
`.venv\Scripts\activate` instead. For the chemistry library alone, install
with `python -m pip install -e .`.

The old launcher remains available:

```bash
python tools/live_renderer.py
```

For an ASGI server or Render deployment:

```bash
python -m pip install -r requirements-render.txt
uvicorn pyPept.web.app:app --host 0.0.0.0 --port 8732
```

Render can also run `python tools/live_renderer.py`, which reads its `PORT`
environment variable. `/health` reports whether the server responds.

## Build a peptide

```python
from rdkit import Chem
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

sequence = Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')
molecule = Molecule(sequence).get_molecule(fmt='ROMol')
print(Chem.MolToSmiles(molecule))
```

`-` always connects the left monomer's R2 to the right monomer's R1.
For standard amino acids, R1 is the backbone nitrogen, R2 the backbone
carboxyl, R3 an additional nitrogen attachment, and R4+ sidechain sites.
Inspect a monomer's actual slots in the library before connecting it.
Some legacy library records use older slot assignments.

| Purpose | Example |
| --- | --- |
| Linear peptide | `A-G-K` |
| Terminal caps | `ac-A-G-K-am` |
| Inline protection | `fmoc-C.trt(4,2)-am` |
| Disulfide | `C.!1(4,4)-A-G-C.!1` |
| Head-to-tail cycle | `!1-A-G-K-A-!1` |
| Sequential branch | `K.[G(4,2).ac(1,2)]-A` |
| Protected branch notation | `K.{G(4,2).ac(1,2)}-A` |

See [the notation guide](docs/notation.md) for branch direction,
attachment slots, and format conversion.

## Convert and compare

```python
from pyPept.smiles import smiles_to_cabiln_core

cabiln, residues = smiles_to_cabiln_core('N[C@@H](C)C(=O)NCC(=O)O')
print(cabiln)
```

Conversion checks the assembled result against the source molecular graph.
Unsupported structures raise an error. When a library match cannot preserve
the input, supported structures can use an inline SMILES fragment instead.
Local synthetic residues retain the surrounding recognized residues when a
verified decomposition is possible. A coarse fallback preserves the component
but cannot identify individual residues inside it; the result reports this.

`pyPept.smiles.convert_smiles` exposes recognition as `complete`, `partial`, or
`unresolved`, with source atom assignments for the identified occurrences.
It also reports when search limits prevent checking further interpretations.
Recognition uses the current monomer library, including newly registered records.
See [the decomposition contract](docs/decomposition.md) for the algorithm and limits.

Unspecified input stereochemistry can be filled from library monomers.
The converter warns when this happens; the web app displays that warning.
Defined stereochemistry and disconnected molecular components must be preserved.
Structure comparison distinguishes an exact match from compatible, incomplete
input stereochemistry. It does not treat different tautomers as exact matches.
A successful conversion does not establish a compound's identity or biological activity.

## Monomer registration

The web app serves a read-only library by default. For a trusted local instance:

```bash
cabiln --enable-registration
```

This enables the Register page and writes to the installed monomer library.
`CABILN_ENABLE_REGISTRATION=1` enables the same behavior for ASGI deployments.
This switch is for a trusted instance; it does not provide user authentication.
Keep public deployments read-only. Back up custom monomers before replacing an
installation or redeploying to an ephemeral filesystem.

Preview detects attachment sites before registration. The server checks that
slots and metadata agree, rejects duplicate symbols, and replaces the SDF
atomically under a file lock. This avoids partial writes and concurrent
registrations with the same symbol.

To keep a custom library outside the installation, copy the distributed SDF
and select it before starting the application or CLI:

```bash
export CABILN_MONOMER_LIBRARY=/path/to/monomers.sdf
cabiln --enable-registration
```

The file must exist. Put a companion `monomers.csv` beside it to retain custom
synonyms. CLI ingestion and web registration use the same selected library.
New records are discovered without a server restart; reopen the palette or
return to its browser window to refresh the tiles. Reaction filters and slot
buttons use the same chemistry detection as assembly. Adding an unsupported
reaction still requires a reaction definition and detection rules.

## Development

```bash
python -m pip install -e '.[dev,web]'
python -m pytest -m 'not distribution'
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
```

The distribution test builds an sdist and wheel, then installs them in a fresh
environment outside the checkout. It needs package-index access. CI checks
Python 3.9, 3.11, and 3.13. Historical repair scripts under `tools/` are not test
entry points; some modify library files when run.

The code layout and review findings are in [docs/architecture.md](docs/architecture.md)
and [docs/cleanup-review.md](docs/cleanup-review.md).

## Attribution

Original pyPept authors: Rodrigo Ochoa, J. B. Brown, and Thomas Fox.
CABILN extensions: Cameron Beeley. See [LICENSE](LICENSE) and
[LICENSE-ORIGINAL](LICENSE-ORIGINAL).

- [pyPept publication](https://doi.org/10.1186/s13321-023-00748-2), Journal of Cheminformatics, 2023.
- [BILN publication](https://doi.org/10.1021/acs.jcim.2c00703), Journal of Chemical Information and Modeling, 2022.
