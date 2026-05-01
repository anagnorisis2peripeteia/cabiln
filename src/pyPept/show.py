"""
pyPept.show — render peptide sequences as 2D skeletal structures.

Tiered renderer, best-available-first:

  ascii      Pure-Python PIL -> block characters.  Always available.
             Warnings: crowded above 10 residues; atom labels invisible.
  terminal   chafa / viu / imgcat — coloured block characters in-terminal.
             Warns if no renderer found, falls back to ascii.
  png        RDKit -> PNG.  Saves to file, optionally opens with OS viewer.
  svg        RDKit -> SVG.  Saves to file.
  browser    SVG in a temp HTML file, opened in default browser.
  mol        MDL MOL block — opens in ChemDraw, Marvin, Avogadro, etc.
  cdxml      ChemDraw XML (coordinates from RDKit 2D layout).

Usage:

    from pyPept.show import show

    show('A-G-K')                     # auto: terminal if chafa present, else ascii
    show('A-G-K', fmt='ascii')        # force ASCII
    show('A-G-K', fmt='png', output='mypep.png')
    show('A-G-K', fmt='browser')      # opens default browser
    show('A-G-K', fmt='mol', output='mypep.mol')   # open in ChemDraw

    # Also accepts an assembled RDKit mol directly:
    from rdkit import Chem
    show(Chem.MolFromSmiles('NCC(=O)O'), fmt='ascii')
"""

__credits__ = ["Cameron Beesley"]
__license__ = "MIT"

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import warnings
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor


# ── Thresholds and constants ────────────────────────────────────────────────

_ASCII_WARN_RESIDUES = 8    # warn about crowding
_ASCII_MAX_RESIDUES  = 20   # suggest switching format

_TERMINAL_RENDERERS = [
    ('chafa', ['chafa', '--size=80x40']),
    ('viu',   ['viu',  '-w', '80']),
    ('imgcat', ['imgcat']),
]


# ── Public entry point ───────────────────────────────────────────────────────

def show(source, fmt='auto', output=None, open_after=True, size=(600, 400)):
    """Render a peptide sequence or molecule.

    :param source: BILN string, pyPept.Sequence, pyPept.Molecule, or RDKit ROMol.
    :param fmt: 'auto' | 'ascii' | 'terminal' | 'png' | 'svg' | 'browser' | 'mol' | 'cdxml'
    :param output: output file path (for png / svg / mol / cdxml).  If None,
                   a temp file is used for file-based formats.
    :param open_after: if True, attempt to open the file with the OS viewer
                       (for png / svg / browser / mol / cdxml).
    :param size: (width, height) in pixels for raster output.
    """
    rdmol, n_residues = _to_rdmol(source)
    if rdmol is None:
        print("[pyPept.show] Could not convert input to a molecule — nothing to draw.")
        return

    rdDepictor.Compute2DCoords(rdmol)

    if fmt == 'auto':
        fmt = _detect_best_fmt()

    _RENDERERS = {
        'ascii':    _render_ascii,
        'terminal': _render_terminal,
        'png':      _render_png,
        'svg':      _render_svg,
        'browser':  _render_browser,
        'mol':      _render_mol,
        'cdxml':    _render_cdxml,
    }

    if fmt not in _RENDERERS:
        print(f"[pyPept.show] Unknown fmt={fmt!r}. "
              f"Choose from: {', '.join(_RENDERERS)}.")
        return

    _RENDERERS[fmt](rdmol, n_residues=n_residues, output=output,
                    open_after=open_after, size=size)


# ── Input normalisation ──────────────────────────────────────────────────────

def _to_rdmol(source):
    """Return (RDKit ROMol, n_residues_estimate).  n_residues may be None."""
    # Already an RDKit mol
    if isinstance(source, Chem.rdchem.Mol):
        return source, None

    # pyPept Molecule wrapper
    try:
        from pyPept.molecule import Molecule as _Mol
        if isinstance(source, _Mol):
            return source.get_molecule(fmt='ROMol'), None
    except ImportError:
        pass

    # pyPept Sequence -> assemble
    try:
        from pyPept.sequence import Sequence as _Seq
        from pyPept.molecule import Molecule as _Mol
        if isinstance(source, _Seq):
            n = source.length()
            mol = _Mol(source).get_molecule(fmt='ROMol')
            return mol, n
    except ImportError:
        pass

    # String — try as BILN first, fall back to SMILES
    if isinstance(source, str):
        # SMILES heuristic
        _SMILES_CHARS = set('()[]=#@+\\/')
        if _SMILES_CHARS & set(source):
            mol = Chem.MolFromSmiles(source)
            return mol, None

        # Try BILN
        try:
            from pyPept.sequence import Sequence as _Seq
            from pyPept.molecule import Molecule as _Mol
            seq = _Seq(source)
            n = seq.length()
            mol = _Mol(seq).get_molecule(fmt='ROMol')
            return mol, n
        except Exception as exc:
            print(f"[pyPept.show] Could not parse {source!r} as BILN: {exc}")
            mol = Chem.MolFromSmiles(source)
            return mol, None

    return None, None


# ── Format detection ─────────────────────────────────────────────────────────

def _detect_best_fmt():
    for name, cmd in _TERMINAL_RENDERERS:
        if shutil.which(cmd[0]):
            return 'terminal'
    return 'ascii'


def _find_terminal_renderer():
    for name, cmd in _TERMINAL_RENDERERS:
        if shutil.which(cmd[0]):
            return name, cmd
    return None, None



# ── Renderers ────────────────────────────────────────────────────────────────

def _render_ascii(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Render as ASCII block art in the terminal."""
    from PIL import Image as _Image

    if n_residues and n_residues >= _ASCII_MAX_RESIDUES:
        warnings.warn(
            f"[pyPept.show/ascii] Sequence has ~{n_residues} residues — ASCII output "
            "will be extremely crowded and mostly unreadable.\n"
            "  Recommended: show(..., fmt='png') or show(..., fmt='browser')",
            UserWarning, stacklevel=3,
        )
    elif n_residues and n_residues >= _ASCII_WARN_RESIDUES:
        warnings.warn(
            f"[pyPept.show/ascii] Sequence has ~{n_residues} residues — ASCII output "
            "shows bond topology but atom labels will be lost.\n"
            "  For readable output try: show(..., fmt='png') or show(..., fmt='browser')",
            UserWarning, stacklevel=3,
        )

    img = Draw.MolToImage(rdmol, size=size)
    img = img.convert('L')

    cols = 72
    aspect = 0.45  # terminal chars are taller than wide
    w, h = img.size
    new_h = int(cols * (h / w) * aspect)
    img = img.resize((cols, new_h))

    chars = '@%#*+=-:. '
    lines = []
    for y in range(img.height):
        row = ''
        for x in range(img.width):
            p = img.getpixel((x, y))
            row += chars[int(p / 256 * len(chars))]
        lines.append(row)

    print()
    print('\n'.join(lines))
    print()
    print("  [ascii renderer — atom labels not visible; "
          "use fmt='png' or fmt='browser' for full quality]")


def _render_terminal(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Render via chafa / viu / imgcat — coloured block characters."""
    name, cmd = _find_terminal_renderer()
    if name is None:
        warnings.warn(
            "[pyPept.show/terminal] No terminal image renderer found.\n"
            "  Install one of:\n"
            "    Linux/Mac:  brew install chafa   or   apt install chafa\n"
            "    Any:        cargo install viu\n"
            "    iTerm2:     imgcat (built-in)\n"
            "  Falling back to ASCII.",
            UserWarning, stacklevel=3,
        )
        _render_ascii(rdmol, n_residues=n_residues, output=output,
                      open_after=open_after, size=size)
        return

    if n_residues and n_residues >= _ASCII_MAX_RESIDUES:
        warnings.warn(
            f"[pyPept.show/terminal] Sequence has ~{n_residues} residues — "
            "terminal rendering will be cramped.\n"
            "  Recommended: show(..., fmt='png') or show(..., fmt='browser')",
            UserWarning, stacklevel=3,
        )

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tf:
        tmp = tf.name

    try:
        _write_png(rdmol, tmp, size=size)
        subprocess.run(cmd + [tmp], check=False)
    finally:
        os.unlink(tmp)

    print(f"\n  [terminal renderer: {name} — "
          "use fmt='png' or fmt='browser' for full quality]")


def _render_png(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Save PNG and optionally open with OS viewer."""
    path = output or tempfile.mktemp(suffix='.png')
    _write_png(rdmol, path, size=size)
    print(f"[pyPept.show] PNG saved -> {path}")
    if open_after:
        _os_open(path)


def _render_svg(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Save SVG and optionally open."""
    path = output or tempfile.mktemp(suffix='.svg')
    w, h = size
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.DrawMolecule(rdmol)
    d.FinishDrawing()
    Path(path).write_text(d.GetDrawingText(), encoding='utf-8')
    print(f"[pyPept.show] SVG saved -> {path}")
    if open_after:
        _os_open(path)


def _render_browser(rdmol, *, n_residues=None, output=None, open_after=True, size=(800, 600)):
    """Wrap SVG in HTML and open in default browser."""
    import webbrowser
    w, h = size
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.DrawMolecule(rdmol)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    smi = Chem.MolToSmiles(rdmol)
    html = textwrap.dedent(f"""\
        <!DOCTYPE html><html><head>
        <meta charset="utf-8">
        <title>pyPept structure</title>
        <style>body{{background:#fff;display:flex;flex-direction:column;
                     align-items:center;font-family:monospace;padding:2em}}</style>
        </head><body>
        {svg}
        <p style="color:#888;font-size:0.8em">SMILES: {smi}</p>
        </body></html>
    """)

    path = output or tempfile.mktemp(suffix='.html')
    Path(path).write_text(html, encoding='utf-8')
    print(f"[pyPept.show] HTML saved -> {path}")
    webbrowser.open(f'file://{path}')


def _render_mol(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Save MDL MOL file. Open in ChemDraw, Marvin Sketch, or Avogadro."""
    path = output or tempfile.mktemp(suffix='.mol')
    molblock = Chem.MolToMolBlock(rdmol, includeStereo=True)
    Path(path).write_text(molblock, encoding='utf-8')
    print(f"[pyPept.show] MOL file saved -> {path}")
    print("  Open with: ChemDraw, Marvin Sketch, Avogadro, or any viewer that reads MDL MOL.")


def _render_cdxml(rdmol, *, n_residues=None, output=None, open_after=True, size=(600, 400)):
    """Save ChemDraw XML (.cdxml) from RDKit 2D coordinates."""
    path = output or tempfile.mktemp(suffix='.cdxml')

    conf = rdmol.GetConformer()
    atoms = rdmol.GetAtoms()
    bonds = rdmol.GetBonds()

    # Build atom id map (RDKit idx -> cdxml id, 1-based)
    aid = {a.GetIdx(): i + 1 for i, a in enumerate(rdmol.GetAtoms())}

    # Scale: RDKit coords are in Angstrom (~1.5 Å bond), cdxml uses ~30 units
    scale = 30.0 / 1.52

    atom_lines = []
    for atom in rdmol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        x = pos.x * scale
        y = -pos.y * scale  # cdxml y-axis is flipped
        sym = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        charge_attr = f' Charge="{charge}"' if charge else ''
        atom_lines.append(
            f'  <n id="{aid[idx]}" p="{x:.2f} {y:.2f}" '
            f'Element="{sym}"{charge_attr}/>'
        )

    bond_order_map = {
        Chem.rdchem.BondType.SINGLE: '1',
        Chem.rdchem.BondType.DOUBLE: '2',
        Chem.rdchem.BondType.TRIPLE: '3',
        Chem.rdchem.BondType.AROMATIC: '1.5',
    }
    bond_lines = []
    for bond in rdmol.GetBonds():
        b = bond_order_map.get(bond.GetBondType(), '1')
        bond_lines.append(
            f'  <b B="{aid[bond.GetBeginAtomIdx()]}" '
            f'E="{aid[bond.GetEndAtomIdx()]}" Order="{b}"/>'
        )

    cdxml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE CDXML SYSTEM "http://www.camsoft.com/xml/cdxml.dtd">
        <CDXML>
         <page>
          <fragment>
        """) + '\n'.join(atom_lines) + '\n' + '\n'.join(bond_lines) + textwrap.dedent("""
          </fragment>
         </page>
        </CDXML>
        """)

    Path(path).write_text(cdxml, encoding='utf-8')
    print(f"[pyPept.show] CDXML saved -> {path}")
    print("  Open with: ChemDraw (File > Open, or double-click if .cdxml is registered).")


# ── Helpers ──────────────────────────────────────────────────────────────────

def _write_png(rdmol, path, size=(600, 400)):
    w, h = size
    try:
        d = rdMolDraw2D.MolDraw2DCairo(w, h)
        d.DrawMolecule(rdmol)
        d.FinishDrawing()
        Path(path).write_bytes(d.GetDrawingText())
    except Exception:
        img = Draw.MolToImage(rdmol, size=size)
        img.save(path)


def _os_open(path):
    """Open a file with the OS default application."""
    try:
        if sys.platform == 'win32':
            os.startfile(path)
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', path])
        else:
            subprocess.Popen(['xdg-open', path])
    except Exception as exc:
        print(f"  (Could not auto-open: {exc})")
