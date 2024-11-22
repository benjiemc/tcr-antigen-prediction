"""Crop TCR-pMHC to binding interface (TCR Variable domain and MHC antigen binding domain)."""

import argparse
import logging
import sys

from Bio.PDB import PDBIO, PDBParser

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.structure import NonHetSelect, crop_structure, extract_chains

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to the pdb')
parser.add_argument('--output', '-o', required=True, help='path to the output pdb')
parser.add_argument('--tcr-chains', nargs='*', help='pdb chains for the tcr structure')
parser.add_argument('--mhc-chains', nargs='*', help='pdb chains for the mhc structure')
parser.add_argument('--antigen-chain', nargs='?', help='pdb chain for the antigen')
parser.add_argument('--remove-het-atoms', action='store_true', help='remove hetero atoms from structure')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    collated_chains = [chain for chain in (args.tcr_chains + args.mhc_chains + [args.antigen_chain]) if chain]

    logger.info('Loading structure')
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('', args.structure)
    structure = extract_chains(structure, collated_chains)

    logger.info('Cropping structure')
    cropped_structure = crop_structure(
        structure, args.tcr_chains, args.mhc_chains, 'MH1' if len(args.mhc_chains) == 1 else 'MH2'
    )

    logger.info('Saving structure')
    io = PDBIO()
    io.set_structure(cropped_structure)

    if args.remove_het_atoms:
        logger.info('Removing hetero atoms')
        logger.info('Saving structure')
        io.save(args.output, NonHetSelect())

    else:
        logger.info('Saving structure')
        io.save(args.output)


if __name__ == '__main__':
    main()
