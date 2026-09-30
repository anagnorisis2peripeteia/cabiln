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

## Use the builder

Open Library and choose **Use** beside a monomer to start a peptide. Select a
residue chip, choose another monomer, select an attachment site on each side,
and press **Connect**. The selected sites light up in their previews. Existing
tile insertion, right-click selection, and selection through the drawing remain
available. The same controls work with newly registered monomers.

Sequence edits, connections, and conversions support **Undo** and **Redo**.
Use Ctrl/Cmd+Z and Ctrl/Cmd+Shift+Z from the sequence input or builder controls.
Changing the input format retains that format's draft and conversion warnings.
The drawing stays visible during updates, with an explicit previous-drawing
notice; selection and exports resume after the current input is validated.

Drafts are saved in this browser. After reloading, choose **Restore saved draft**
to recover the sequence, format drafts, and reference text. Starting new work
resolves the recovery offer and saves the new draft. Browser storage must be
available for recovery; Undo still works when storage is unavailable.

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
Pass `notation="bracket"` to `convert_smiles` for bracket output; its source atom
assignments refer directly to occurrences in that output. Editing, formatting
and highlighting share explicit occurrence identities and numbered connections.
The [independent benchmark](docs/decomposition-benchmark.md) checks molecular
identity, residue boundaries and attachment sites against external references.

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

This enables the Register page for loopback clients. Use an external library
for durable custom data. `CABILN_ENABLE_REGISTRATION=1` enables the same route
for ASGI deployments. For remote administration, set a secret
`CABILN_REGISTRATION_TOKEN` of at least 32 characters; the browser prompts for
HTTP Basic credentials (username `admin`, or `CABILN_REGISTRATION_USER`). Serve
remote administration over HTTPS. The public builder stays read-only by default.
Production administration also requires an external library and backup directory.

Preview detects attachment sites before registration. The server checks that
slots and metadata agree, rejects duplicate symbols, and replaces the SDF
atomically under a file lock. This avoids partial writes and concurrent
registrations with the same symbol.

To keep a custom library outside the installation, copy the distributed SDF
and select it before starting the application or CLI:

```bash
export CABILN_MONOMER_LIBRARY=/path/to/monomers.sdf
export CABILN_LIBRARY_BACKUP_DIR=/path/to/backups
cabiln --enable-registration
```

The file must exist. Put a companion `monomers.csv` beside it to retain custom
synonyms. CLI ingestion and web registration use the same selected library.
New records are discovered without a server restart; reopen the palette or
return to its browser window to refresh the tiles. Reaction filters and slot
buttons use the same chemistry detection as assembly. Adding an unsupported
reaction still requires a reaction definition and detection rules.

Each registration snapshots the existing SDF and companion aliases when the
backup directory is configured. If that backup fails, registration does not
write. [Deployment and recovery](docs/deployment.md) describes offline restore.
The library's [quality baseline](docs/library-quality.md) identifies known
compatibility exceptions without changing stored structures or attachment slots.

## Saved work and production

Save project downloads editable source, notation drafts, original reference,
recognition details, and the library/rule binding. Open project checks that its
sources and chemistry still agree before replacing current work. Adding an
unrelated monomer can remain compatible; changed selected definitions or rules
require re-verification. Browser drafts remain local, with a clear-draft control
in Help. Server processing and local storage are explained there.

Production uses a bounded worker process with a deadline, cancellation cleanup,
memory limit on Linux, and explicit overload responses. `/health` reports
liveness; `/ready` checks startup data and worker availability. The container
profile freezes Python and runtime dependencies and disables registration.
See [deployment](docs/deployment.md) and the current
[launch evidence](docs/launch-validation.md) before deploying.

## Development

```bash
python -m pip install -e '.[dev,web]'
python -m pytest -m 'not distribution and not fuzz'
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
```

The distribution test builds an sdist and wheel, then installs them in a fresh
environment outside the checkout. It needs package-index access. CI checks
Python 3.9, 3.11, and 3.13. Historical repair scripts under `tools/` are not test
entry points; some modify library files when run.

Generated tests use Hypothesis for notation, chemistry and library changes, and
fast-check for UI histories and browser interactions. Run the bounded campaign:

```bash
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/python python -m pytest -m fuzz \
  --hypothesis-show-statistics --timeout=300 --timeout-method=thread
npm ci --prefix tests/browser
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/frontend npm run fuzz:frontend --prefix tests/browser
# Install Chromium first: cd tests/browser && npx playwright install chromium
CABILN_FUZZ_ARTIFACTS=/tmp/cabiln-fuzz/browser npm run fuzz:browser --prefix tests/browser
```

Set `CABILN_FUZZ_PROFILE=deep` for longer campaigns; use `--timeout=1200` for the
Python run. CI runs bounded campaigns on pushes and pull requests, and the deep
profile nightly. Browser servers and library changes use temporary local copies.
These campaigns sample supported feature combinations; they do not exhaust the
possible peptides or editing histories.

Failures retain reduced inputs, seeds, dependency versions and replay information.
Python artifacts include source/library hashes, Hypothesis's example database
and reproduction decorator. Replay with the recorded versions and
`--hypothesis-seed=<seed>`, or apply the reported reproduction decorator to the
owning test. Select the recorded fast-check property with `CABILN_FUZZ_CASE`, then
use `CABILN_FUZZ_SEED`, `CABILN_FUZZ_PATH` and, for command histories,
`CABILN_FUZZ_REPLAY_PATH` to replay its failure. `CABILN_FUZZ_SCHEDULE` accepts the
saved task order for exact replay of a frontend scheduling failure.
Observation counters include shrinking and repeated examples; they do not count
unique peptides. A watchdog termination preserves
the latest eight Python observations in `active.json`, alongside the captured CI log.
The [initial campaign record](tools/benchmarks/results/fuzz-20260930.json) documents
generated cases, independent fault controls and coverage limits.

The code layout and review findings are in [docs/architecture.md](docs/architecture.md)
and [docs/cleanup-review.md](docs/cleanup-review.md). The latest editing, browser,
and performance evidence is in
[docs/ux-performance-validation.md](docs/ux-performance-validation.md).

## Attribution

Original pyPept authors: Rodrigo Ochoa, J. B. Brown, and Thomas Fox.
CABILN extensions: Cameron Beeley. See [LICENSE](LICENSE) and
[LICENSE-ORIGINAL](LICENSE-ORIGINAL).

- [pyPept publication](https://doi.org/10.1186/s13321-023-00748-2), Journal of Cheminformatics, 2023.
- [BILN publication](https://doi.org/10.1021/acs.jcim.2c00703), Journal of Chemical Information and Modeling, 2022.
