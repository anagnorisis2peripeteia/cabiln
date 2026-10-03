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
from pathlib import Path

from rdkit import Chem

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from pyPept.monomer_store import library_path, monomer_record, register_molecule
from pyPept.interfaces.monomer_pipeline import pre_activate, ActivationError


def _default_sdf_path():
    """Return the configured default monomer library."""
    return library_path()


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


def register_monomer(smiles, symbol, name=None, m_type='aa', m_subtype='modified',
                     sdf_path=None, backbone_indices=None):
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
        configured default library when ``None``.
    :param backbone_indices: {1: n_idx, 2: c_idx} to force specific backbone
        atoms for β/γ orientation variants.  Passed to
        :func:`~pyPept.interfaces.monomer_pipeline.pre_activate` unchanged.
    :returns: :class:`~pyPept.interfaces.monomer_pipeline.ActivationResult`.
    :raises ActivationError: if :func:`pre_activate` cannot process the SMILES.
    :raises ValueError: if the symbol or an alias already exists in the library.
    """
    if sdf_path is None:
        sdf_path = _default_sdf_path()
    sdf_path = Path(sdf_path)
    if name is None:
        name = symbol

    result = pre_activate(smiles, backbone_indices=backbone_indices)

    mol = monomer_record(
        result.chuckles, symbol, result.leaving, result.chem_types,
        name=name, m_type=m_type, m_subtype=m_subtype,
        activation_policy=result.policy,
    )

    register_molecule(mol, sdf_path=sdf_path)
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
             "Defaults to CABILN_MONOMER_LIBRARY or the installed monomers.sdf.")
    parser.add_argument(
        '--backbone', nargs=2, type=int, metavar=('R1_ATOM', 'R2_ATOM'),
        help='Select a backbone using zero-based atom indices in the input SMILES.')

    args = parser.parse_args()

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

    try:
        result = register_monomer(
            smiles,
            symbol=args.symbol,
            name=args.name,
            m_type=args.m_type,
            m_subtype=args.m_subtype,
            sdf_path=sdf_path,
            backbone_indices=dict(zip((1, 2), args.backbone)) if args.backbone else None,
        )
    except ActivationError as exc:
        print(f"ERROR: Pre-activation failed: {exc}", file=sys.stderr)
        sys.exit(1)
    except (ValueError, OSError) as exc:
        print(f"ERROR: Registration failed: {exc}", file=sys.stderr)
        sys.exit(1)

    target = sdf_path if sdf_path else _default_sdf_path()
    print(f"Registered '{args.symbol}' -> {target}")
    print(f"  CHUCKLES  : {result.chuckles}")
    print(f"  Leaving   : {result.leaving}")
    if result.chem_types:
        print(f"  Chem types: {result.chem_types}")


if __name__ == "__main__":
    main()
