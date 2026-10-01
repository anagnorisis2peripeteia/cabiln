"""
Command line interface to interact with pyPept

From the publication: 
'pyPept: a python library to generate atomistic 2D and 3D representations of peptides'
Journal of Cheminformatics, 2023
"""


__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"
__version__ = "1.0"


import argparse
import logging
import sys
import warnings

from rdkit import Chem
from rdkit.Chem import Draw

from pyPept.inputs import read_input
from pyPept.pdb_names import correct_pdb_atoms
from pyPept.conformer import Conformer, ConformerConstants, SecStructPredictor


_def_image_size = (1200, 1200)

def get_inputs_parser():
    """Constructs parser for inputs."""

    parser = argparse.ArgumentParser(add_help=False)

    input_type_group = parser.add_mutually_exclusive_group(required=True)
    input_type_group.add_argument(
        '--biln', type=str, metavar='string',
        required=False,
        help="BILN string with the peptide to analyze.")
    input_type_group.add_argument(
        '--helm', type=str, metavar='string',
        required=False,
        help="HELM string with the peptide to analyze.")
    input_type_group.add_argument(
        '--fasta', type=str, metavar='string',
        required=False,
        help="FASTA string with the peptide to analyze.\
              Only natural amino acids are allowed.")

    additional_type_group = parser.add_argument_group('Additional options')
    additional_type_group.add_argument(
        '--depiction', type=str, metavar='text',
        required=False, default='local',
        help="""Method to generate the 2D image. 
            Two options are supported: 'local' (default) or 'rdkit'.
            If the default rendering does not result in a horizontal
            alignment easy for viewing, choose the RDKit-based rendering
            option (such as for the molecule PEPTIDEPEPTIDEPEPTIDE).
            """
            )
    additional_type_group.add_argument(
        '--prefix', type=str, metavar='text',
        required=False, default='peptide',
        help="Name used in the output files. The default is 'peptide'.")
    additional_type_group.add_argument(
        '--secstruct', type=str, metavar='text',
        required=False, default=None,
        help="Use the given secondary structure. " + \
             "Otherwise, the secondary structure is predicted and used.")
    additional_type_group.add_argument(
        '--sdf2D', action="store_true",
        help="Generate a 2D SDF file of the peptide.")
    additional_type_group.add_argument(
        '--noconf', action="store_true",
        help="Do not generate a conformer for the peptide.")
    additional_type_group.add_argument(
        '--imagesize', type=int, metavar='dim', nargs=2,
        required=False, default=_def_image_size,
        help="Image size for 2D depiction, default %s." % (_def_image_size,))

    log_options = parser.add_argument_group('Logging options')
    log_options.add_argument(
        '--logfile', type=str, metavar='filename',
        required=False, default=None,
        help="Output messages to given logfile, default is stderr.")
    log_options.add_argument(
        "-v", "--verbose", action="store_true",
        help="Increase output verbosity")

    return parser

def main():
    useParser = argparse.ArgumentParser(
        description="""Generate atomistic 2D and 3D representations of peptides
        from given monomer sequences.""",
        parents=(get_inputs_parser(),
                 ))
    args = useParser.parse_args()

    log_stream = logging.StreamHandler(
        open(args.logfile, "w", encoding='utf-8') 
        if args.logfile else sys.stderr)
    log_stream.setFormatter(logging.Formatter(
        "%(asctime)-10s %(levelname)s:%(message)s",
        datefmt="%H:%M:%S"))
    logger = logging.Logger(name="pyPept")
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logger.addHandler(log_stream)
    logger.debug("Invocation arguments: %s", args)

    outFileList = list()

    selected = next(
        ((kind, getattr(args, kind)) for kind in ("biln", "helm", "fasta")
         if getattr(args, kind)),
        None,
    )
    if selected is None:
        logger.error(
            "An input should be provided using BILN, HELM or FASTA format.")
        sys.exit(1)
    kind, source = selected
    # The historical --biln flag accepts current CABILN, as before.
    parsed = read_input(
        source, input_format="cabiln" if kind == "biln" else kind,
        warning_sink=warnings.warn,
    )
    biln, seq = parsed.source, parsed.sequence
    label = "BILN" if kind == "biln" else f"{kind.upper()}->BILN"
    logger.info("1. Processing the %s sequence %s", label, biln)

    if not args.noconf:
        parsed.sequence = seq = correct_pdb_atoms(seq)

    logger.info("2. Creating the RDKit object")
    if args.depiction in ["local", "rdkit"]:
        romol = parsed.assemble(depiction=args.depiction)
        mol = parsed.assembly
    else:
        logger and logger.error(
            "Please select a depiction mode from (local, rdkit)")
        exit(2)
    print(f"The SMILES of the peptide is: {Chem.MolToSmiles(romol)}")
    if args.sdf2D:
        mol.write_molecule(fmt='SDF', out_file=f'{args.prefix}.sdf')

    # Render peptide to an image file.
    outImage = f'{args.prefix}.png'
    Draw.MolToFile(romol, outImage, size=args.imagesize)
    outFileList.append(outImage)

    # Create the peptide conformer with correct atom names and SS
    if not args.noconf:
        logger.info(
            "3. Predicting the peptide conformer")
        fasta = Conformer.get_peptide(biln)
        if args.secstruct is None:
            ss_input = SecStructPredictor.predict_active_ss(fasta)
        else:
            invalid = [v for v in args.secstruct
                       if v not in ConformerConstants.expected_ss_symbols]
            if len(invalid) > 0:
                raise RuntimeError(
                    "{} invalid secondary structure symbols (not {}): {} .".format(
                        len(invalid),
                        " ".join(ConformerConstants.expected_ss_symbols),
                        " ".join(invalid)))
            if len(args.secstruct) == len(fasta):
                ss_input = args.secstruct
            else:
                logger.error(
                    "Check the length of the input secondary structure. " + \
                    "It should be the same than the peptide main chain " + \
                    "(without capping groups and extensions). " + \
                    f"Peptide input length: {len(fasta)}, " + \
                    f"Secondary Structure input length: {len(args.secstruct)}."
                    )
                exit(3)
            print(f"The provided Secondary Structure is: {ss_input}")
        romol = Conformer.generate_conformer(romol, ss_input, 
                    generate_pdb=True, output_name=args.prefix)
        outFileList.append(f'{args.prefix}.pdb')
    
    for f in outFileList:
        logger.info(f"File generated: {f}.")

if __name__ == "__main__":
    main()
