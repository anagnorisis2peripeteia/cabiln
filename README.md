# CABILN / pyPept

CABILN is a peptide notation and interactive builder based on
[Boehringer Ingelheim's pyPept](https://github.com/Boehringer-Ingelheim/pyPept).
Its monomer library supplies the builder's vocabulary: registration detects
attachment sites and chemistry, creates a library tile, and makes the monomer
available for construction and recognition. New monomers with supported
attachment chemistry need no code specific to their names.

The builder supports backbones, caps, branches, crosslinks, monomer replacement,
structure comparison, and PNG/MOL export. Monomer structures and reaction
definitions determine the assembled RDKit molecule.

**[Open the web app](https://cabiln.onrender.com/)** · [Documentation](docs/README.md)

## Install and run

Use Python 3.9 or newer. From a checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[web]'
cabiln
```

Open `http://127.0.0.1:8732`. On Windows, activate the environment with
`.venv\Scripts\activate`. For the chemistry library alone, install with
`python -m pip install -e .`.

`cabiln --host` and `--port` set the listening address; the launcher also reads
`HOST` and `PORT`. `python tools/live_renderer.py` remains supported.
See [deployment](docs/deployment.md) for hosting, process limits and recovery.

## Use the builder

Open **Tutorial** in the top toolbar for lessons on building a peptide and
swapping a building block in Retatrutide. It opens a separate practice tab.
Follow the arrow or highlighted button for each action. The lessons mark the exact
tile, connection point or control, then confirm success and point to **Next**.
Open **Why this step?** for the explanation; you do not need to recognise the
chemical structure. Brief animations settle into steady cues, and reduced-motion
settings keep the cues still.
Practice edits do not read or replace your saved browser draft; use **Save
project** if you want to keep them. Browse or search for your chosen amino acid;
the guide follows the resulting product. **Show control** locates the next control
or offers to reopen a closed panel. **Back** reviews instructions without changing
the molecule; **Undo** reverses an edit, and **Restart** reloads the starting
peptide. Closing the guide keeps the tab in practice mode.
On desktop, the guide sits beside the workspace. On smaller screens, **Show
control** brings the target into the space above the guide. After a preview is
ready, **Show preview** brings that drawing into view.

Open **Library** and choose **Use** beside a monomer to start a peptide.
Open **Build**, select a residue tile or click/tap an atom label or bond in the
drawing, choose another monomer, then select a free attachment site on each side.
The selected sites light up in the previews. **Preview** shows the proposed
structure, reaction and CABILN before **Connect** applies the change. You can also
select two existing residues to
connect them, or insert a monomer between neighbouring backbone residues.
The product preview marks the changed block in blue and the new connections in
orange. Open **Reaction and notation** for the reaction details and proposed
source. When only one compatible site pair remains, **Use … R… → … R…** offers
to select it; you still review and apply the connection yourself.

Site numbers stay fixed as you build. The labels and connection checks use the
selected residue's current structure: after N-alkylation, a primary amine's
remaining slot becomes a secondary amine; after N-acylation, it becomes an amide
nitrogen. A free slot is unavailable if another connection consumed its reactive
group. Undo restores the earlier structure and its available chemistry.
Sulfonamide, urea and carbamate nitrogens have separate labels and reaction
rules. Used sulfur sites describe the resulting thioether, thioester or
disulfide. Phosphates expose one port per remaining OH; each stays available
until connected.
These checks describe supported structural transformations, not reaction
conditions or expected yield.

**Filter** in Library opens Build when needed. Select a residue in the drawing or
its tile to match its free attachment sites, then choose a numbered site to narrow
the results. The library states which residue and site it is matching. Turn
Filter off to show all monomers again. Search works within the filtered list.
Category and collection menus narrow the list further. Star a monomer to keep
it in **Favourites**, or choose **Recently used** for your last twelve choices.
These preferences stay in this browser; practice tabs keep separate, temporary
choices. **Reset** clears search and filters without changing the peptide.
**Site details** shows each monomer's attachment codes and leaving groups.

Choose **Swap monomer** in Build to replace a selected residue. The library
shows candidates with compatible sites for its existing connections; choose
**Select** beside a replacement. Review the mapping from each occupied R-group
to a distinct replacement site in **Keep all connections**; different site numbers
are supported. **Change** on the replacement card keeps the original residue
selected while you choose again. **Preview swap** assembles the proposed product
before **Apply swap** becomes available. Changing the replacement or mapping requires
a new preview. Neighbours, branches and ring closures stay connected.
The preview marks the replacement in blue and the retained connections in orange.
**Focus change** centres and enlarges that region; **Reset view** shows the whole peptide.
Site renumbering can change the notation's layout. Chemistry filters
narrow the choices; full assembly can still reject a candidate.
**Learn Swap** opens the Retatrutide lesson directly from Build.

Hover, tap or keyboard-activate residue or branch tiles to highlight their atoms.
Selected residue buttons also expose their selection state to assistive tools. Library tiles
support insertion and right-click selection. Drawing, reference and product
preview canvases have zoom and reset controls, keyboard panning, and touch
pan/pinch. Escape closes the active panel and returns focus to its control.
The drawing starts on a charcoal canvas; **Dark** switches it to white.
**PNG** and **MOL** exports sit beside the drawing controls.

**Undo** and **Redo** cover sequence edits, connections, swaps and conversions.
Use Ctrl/Cmd+Z and Ctrl/Cmd+Shift+Z in the sequence input or builder controls.
Changing input format retains that format's draft and conversion warnings.
During updates, the previous drawing stays visible and is marked as outdated;
selection and exports resume when the current input has a valid drawing.
Unknown monomer errors offer **Select problem** to locate the name in your
sequence and explain how to correct it. Failed connection checks keep your
selections so you can retry.

Browser drafts retain the sequence, format drafts and reference. After reloading,
choose **Restore saved draft** to recover them. **Save project** downloads a
portable file; **Open project** checks its saved definitions and chemistry before
replacing current work. Help groups instructions by task, including storage and
**Clear saved browser draft**. **Undo clear** recovers that copy until the next edit.
Undo remains available when browser storage is unavailable; use **Save project**
to keep the session outside the browser.
If the server cannot finish the project check, Retry repeats the save without
discarding your current work.

## Build a peptide in Python

```python
from rdkit import Chem
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

sequence = Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')
molecule = Molecule(sequence).get_molecule(fmt='ROMol')
print(Chem.MolToSmiles(molecule))
```

`-` connects the left monomer's R2 to the right monomer's R1. Standard amino acids
use R1 for the backbone nitrogen, R2 for the backbone carboxyl, R3 for an
additional nitrogen attachment, and R4+ for sidechains. Other definitions can
use different numbering; inspect their actual sites in the library.

| Purpose | Example |
| --- | --- |
| Linear peptide | `A-G-K` |
| Terminal caps | `ac-A-G-K-am` |
| Inline protection | `fmoc-C.trt(4,2)-am` |
| Disulfide | `C.!1(4,4)-A-G-C.!1` |
| Head-to-tail cycle | `!1-A-G-K-A-!1` |
| Sequential branch | `K.[G(4,2).ac(1,2)]-A` |
| Protected branch | `K.{G(4,2).ac(1,2)}-A` |

The [notation guide](docs/notation.md) explains nested and sibling branches,
attachment slots, and bracket/percent conversion. **Canonical** formatting
produces a stable spelling for a resolved monomer graph under a specific library
and convention; it preserves the chosen decomposition.

## Convert and compare

The input selector accepts CABILN, SMILES, BILN and HELM. Verify compares the
current peptide with reference text or an uploaded MOL/SDF structure.

```python
from pyPept.smiles import convert_smiles

result = convert_smiles('N[C@@H](C)C(=O)NCC(=O)O')
print(result.cabiln)
print(result.recognition_status)
```

The converter uses the selected monomer library and checks its assembled result
against the source structure. It reports recognition as `complete`, `partial`, or
`unresolved`, with source atom assignments and any search limits reached.
Supported unknown regions can use local inline SMILES while retaining nearby
recognized residues. An opaque fallback preserves a component without claiming
residue boundaries inside it. Unsupported structures raise an error.

Pass `notation="bracket"` for bracket output. Returned assignments identify the
occurrences in that output. The older `smiles_to_cabiln_core` interface remains
available. See [recognition](docs/decomposition.md) for the algorithm and limits,
and the [benchmark](docs/decomposition-benchmark.md) for independently specified
molecules, residue partitions and attachment sites.

Conversion preserves defined stereochemistry and disconnected components.
Library monomers can supply stereochemistry omitted by the input; the result
and web app report this. Verify distinguishes an exact match from compatibility
with incomplete input stereochemistry. Different tautomers are not exact matches.

## Register monomers

The public app's library is read-only. Enable registration on a trusted local
instance with:

```bash
cabiln --enable-registration
```

Preview detects attachment sites before registration. Review its numbered markers
beside the naming form; **Generated notation** contains the detected CHUCKLES.
If several main-chain or cap orientations are possible, choose a numbered
drawing first. Registration stays disabled until you choose; editing the input
requires a new preview. Equivalent SMILES produce the same automatic R numbering.
Empty input and failed previews explain how to continue without losing entries.
The server validates the
slots and metadata, rejects duplicate symbols, and replaces the SDF atomically
under a file lock. Loading a library does not run detection or change its numbers.
All 1,128 bundled definitions have been reprocessed under the current policy;
older libraries and notation can be updated with the [migration tool](docs/monomer-migration.md).
See the
[attachment contract](docs/architecture.md#attachment-detection-and-registration)
for numbering, explicit overrides and bulk re-imports.

For durable custom data, copy the distributed `monomers.sdf` outside the
installation and select it before starting the application or CLI:

```bash
export CABILN_MONOMER_LIBRARY=/path/to/monomers.sdf
export CABILN_LIBRARY_BACKUP_DIR=/path/to/backups
cabiln --enable-registration
```

The SDF must already exist. Put a companion `monomers.csv` beside it to retain
custom synonyms. CLI ingestion and web registration use the same selected
library. Reopen the palette or return to its browser window to discover new
records without restarting the server. Reaction filters and attachment buttons
use the same chemistry detection as assembly. Supporting a new reaction can
require additional reaction definitions and detection rules.

A configured backup directory receives a snapshot before each registration;
a failed backup prevents the write. Remote administration requires authentication
and HTTPS; production also requires an external library and backup directory.
See [administration and restore](docs/deployment.md#administrative-ingestion-and-backups).

The [library quality baseline](docs/library-quality.md) records compatibility
exceptions and identity concerns for specific definitions. Successful assembly
alone does not verify a library record's chemical identity.

## Development and checks

```bash
python -m pip install -e '.[dev,web]'
python -m pytest -m 'not distribution and not fuzz'
node --test tests/test_frontend.js
python -m pytest tests/test_distribution.py
```

The distribution check builds an sdist and wheel and installs them outside the
checkout; it needs package-index access. CI tests Python 3.9, 3.11 and 3.13.
The browser suite uses Node 20 or newer; CI uses Node 22 and tests Chromium,
Firefox and WebKit. Automated accessibility scans accompany the interaction tests.

See [browser checks](tests/browser/README.md) for Playwright setup and
[generated tests and replay](tools/README.md#generated-tests) for Hypothesis and
fast-check. Generated campaigns sample feature combinations and preserve reduced
failures; they do not exhaust possible peptides or editing histories.

[Architecture](docs/architecture.md) maps the current modules.
[Documentation](docs/README.md) links the usage guides and performance records.
See the [CI runs](https://github.com/anagnorisis2peripeteia/pyPept/actions) for
results from a particular release.

## Attribution

Original pyPept authors: Rodrigo Ochoa, J. B. Brown, and Thomas Fox.
CABILN extensions: Cameron Beeley. See [LICENSE](LICENSE) and
[LICENSE-ORIGINAL](LICENSE-ORIGINAL).

- [pyPept publication](https://doi.org/10.1186/s13321-023-00748-2), Journal of Cheminformatics, 2023.
- [BILN publication](https://doi.org/10.1021/acs.jcim.2c00703), Journal of Chemical Information and Modeling, 2022.
