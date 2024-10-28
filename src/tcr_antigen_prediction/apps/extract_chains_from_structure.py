'''Extract chains from a PDB structue into a new PDB file.'''
import argparse
import logging
import sys

from Bio.PDB import PDBParser, PDBIO

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.structure import extract_chains

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('structure', help='path to pdb file')
parser.add_argument('--output', '-o', required=True, help='path to output pdb file')
parser.add_argument('--chains', nargs='+', required=True, help='chains to extract')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('', args.structure)

    output_structure = extract_chains(structure, args.chains)

    io = PDBIO()
    io.set_structure(output_structure)
    io.save(args.output)


if __name__ == '__main__':
    main()
