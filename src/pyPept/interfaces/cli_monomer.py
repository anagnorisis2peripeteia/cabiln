#!/usr/bin/env python3
"""
CLI to register novel monomers into the pyPept monomer library.

Accepts a plain SMILES or a CABILN sequence (single-residue or small fragment),
runs the pre-activation pipeline, and appends the result to the library SDF.

Usage examples::

    # From plain SMILES
    pyPept-monomer-add --smiles "N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O" \\
        --symbol Lys_Boc --name "Boc-lysine"

    # From CABILN (assemble first, then register)
    pyPept-monomer-add --from-cabiln "K.Boc(4,1)" \\
        --symbol Lys_Boc --name "Boc-lysine"

From publication: pyPept: a python library to generate atomistic 2D and 3D
representations of peptides, Journal of Cheminformatics, 2023.
Updated 2025.
"""

__credits__ = ["J.B. Brown", "Thomas Fox", "Cameron Beesley"]
__license__ = "MIT"

import argparse
import sys
from importlib.resources import files
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from pyPept.interfaces.monomer_pipeline import pre_activate, ActivationError, _LG_COLS


_DEFAULT_SDF_PKG = "pyPept.data"
_DEFAULT_SDF_NAME = "monomers.sdf"


def _default_sdf_path():
    """Return the Path to the installed monomers.sdf."""
    return Path(str(files(_DEFAULT_SDF_PKG).joinpath(_DEFAULT_SDF_NAME)))


def smiles_from_cabiln(cabiln):
    """
    Assemble a CABILN sequence and return the canonical SMILES of the product.

    Designed for single-residue CABILN notation such as ``"K.Boc(4,1)"``.
    The assembled molecule has all R-group dummies replaced by the correct
    leaving-group atoms (OH, H, etc.) — the same path that final peptide
    assembly follows — so the result is a plain, unactivated SMILES suitable
    for passing directly to :func:`pre_activate`.

    :param cabiln: CABILN/BILN string (typically a single residue with caps).
    :returns: canonical SMILES string.
    :raises ValueError: if the CABILN is invalid or assembly yields no molecule.
    """
    seq = Sequence(cabiln)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt='ROMol')
    if romol is None:
        raise ValueError(f"Assembly of CABILN {cabiln!r} produced no molecule.")
    return Chem.MolToSmiles(romol)


def _append_mol_to_sdf(mol, sdf_path):
    """
    Append a single RDKit mol (with properties already set) to an SDF file.

    Writes the V2000 mol block followed by SDF ``>  <prop>`` property tags and
    the ``$$$$`` record terminator, in append mode so existing records are
    preserved.

    :param mol: RDKit Mol with all properties set.
    :param sdf_path: path (str or Path) of the target SDF file.
    """
    mol_block = Chem.MolToMolBlock(mol)
    with open(Path(sdf_path), 'a', encoding='utf-8') as f:
        f.write(mol_block)
        for prop in mol.GetPropNames():
            f.write(f">  <{prop}>\n{mol.GetProp(prop)}\n\n")
        f.write("$$$$\n")


def register_monomer(smiles, symbol, name=None, m_type='aa', m_subtype='modified',
                     sdf_path=None):
    """
    Pre-activate a monomer SMILES and append it to the library SDF.

    This is the programmatic API equivalent of the ``pyPept-monomer-add`` CLI.
    The SMILES is run through :func:`~pyPept.interfaces.monomer_pipeline.pre_activate`
    to derive CHUCKLES, leaving groups, and chem-type annotations, then the
    result is written as a new record in the target SDF.

    :param smiles: plain (non-CHUCKLES) SMILES of the new monomer.
    :param symbol: short token / abbreviation (e.g. ``'Lys_Boc'``).
    :param name: human-readable name; defaults to *symbol* when ``None``.
    :param m_type: monomer type (``'aa'``, ``'cap'``, …).  Default ``'aa'``.
    :param m_subtype: subtype (``'modified'``, ``'natural'``, ``'cap'``, …).
        Default ``'modified'``.
    :param sdf_path: path to the library SDF to append to.  Defaults to the
        installed ``monomers.sdf`` when ``None``.
    :returns: :class:`~pyPept.interfaces.monomer_pipeline.ActivationResult`.
    :raises ActivationError: if :func:`pre_activate` cannot process the SMILES.
    """
    if sdf_path is None:
        sdf_path = _default_sdf_path()
    sdf_path = Path(sdf_path)
    if name is None:
        name = symbol

    result = pre_activate(smiles)

    mol = Chem.MolFromSmiles(result.chuckles)
    if mol is None:
        raise ActivationError(
            f"pre_activate produced unparseable CHUCKLES for {symbol!r}: "
            f"{result.chuckles!r}")

    rdDepictor.Compute2DCoords(mol)
    mol.SetProp('symbol', symbol)
    mol.SetProp('m_abbr', symbol)
    mol.SetProp('m_name', name)
    mol.SetProp('m_type', m_type)
    mol.SetProp('m_subtype', m_subtype)

    leaving = result.leaving
    max_slot = max((s for s, v in leaving.items() if v), default=3)
    lg_vals = [leaving.get(s) for s in range(1, max_slot + 1)]
    mol.SetProp('m_Rgroups', ','.join(str(v) if v else 'None' for v in lg_vals))

    if result.chem_types:
        mol.SetProp('m_chem_types',
                    ','.join(f"{s}:{ct}" for s, ct in sorted(result.chem_types.items())))

    _append_mol_to_sdf(mol, sdf_path)
    return result


def main():
    """Entry point for the ``pyPept-monomer-add`` console script."""

    parser = argparse.ArgumentParser(
        description="Register a novel monomer into the pyPept library SDF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  pyPept-monomer-add --smiles "N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O" --symbol Lys_Boc
  pyPept-monomer-add --from-cabiln "K.Boc(4,1)" --symbol Lys_Boc --name "Boc-lysine"
        """)

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        '--smiles', metavar='SMILES',
        help="Plain SMILES of the new monomer.")
    src.add_argument(
        '--from-cabiln', metavar='CABILN', dest='from_cabiln',
        help="CABILN sequence to assemble (e.g. 'K.Boc(4,1)'); "
             "the resulting plain SMILES is then pre-activated.")

    parser.add_argument(
        '--symbol', required=True, metavar='TOKEN',
        help="Short token / abbreviation for the monomer (e.g. Lys_Boc).")
    parser.add_argument(
        '--name', default=None, metavar='TEXT',
        help="Human-readable name.  Defaults to --symbol.")
    parser.add_argument(
        '--type', default='aa', dest='m_type', metavar='TYPE',
        help="Monomer type ('aa', 'cap', …).  Default: aa.")
    parser.add_argument(
        '--subtype', default='modified', dest='m_subtype', metavar='SUBTYPE',
        help="Monomer subtype ('modified', 'natural', 'cap', …).  Default: modified.")
    parser.add_argument(
        '--sdf', default=None, metavar='PATH',
        help="Path to the library SDF to append to.  "
             "Defaults to the installed monomers.sdf.")

    args = parser.parse_args()

    # --- Resolve input SMILES ---
    if args.from_cabiln:
        try:
            smiles = smiles_from_cabiln(args.from_cabiln)
        except Exception as exc:
            print(f"ERROR: CABILN assembly failed: {exc}", file=sys.stderr)
            sys.exit(1)
        print(f"Assembled SMILES from CABILN: {smiles}")
    else:
        smiles = args.smiles

    sdf_path = Path(args.sdf) if args.sdf else None

    # --- Register ---
    try:
        result = register_monomer(
            smiles,
            symbol=args.symbol,
            name=args.name,
            m_type=args.m_type,
            m_subtype=args.m_subtype,
            sdf_path=sdf_path,
        )
    except ActivationError as exc:
        print(f"ERROR: Pre-activation failed: {exc}", file=sys.stderr)
        sys.exit(1)

    target = sdf_path if sdf_path else _default_sdf_path()
    print(f"Registered '{args.symbol}' -> {target}")
    print(f"  CHUCKLES  : {result.chuckles}")
    print(f"  Leaving   : {result.leaving}")
    if result.chem_types:
        print(f"  Chem types: {result.chem_types}")


if __name__ == "__main__":
    main()
