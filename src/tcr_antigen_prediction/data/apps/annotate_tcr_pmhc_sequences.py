"""Get TCR CDR sequences and peptide sequences for a TCR:pMHC structure. Structure must be IMGT numbered."""

import argparse
import logging
import sys

from Bio.PDB import PDBParser

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3
from tcr_antigen_prediction.data.structure import get_sequence

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to TCR:pMHC pdb structure')
parser.add_argument('--output', '-o', help='output path')

chain_ids_group = parser.add_argument_group('Chain IDs')
chain_ids_group.add_argument('--alpha-chain-id', default=None, help='TCR alpha chain id (Default: None)')
chain_ids_group.add_argument('--beta-chain-id', default=None, help='TCR beta chain id (Default: None)')
chain_ids_group.add_argument('--antigen-chain-id', default=None, help='antigen-chain-id (Default: None)')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('', args.structure)

    header = []
    output = []

    if args.alpha_chain_id:
        for cdr_number, imgt_numbering in enumerate((IMGT_CDR1, IMGT_CDR2, IMGT_CDR3), 1):
            header.append(f'CDR{cdr_number}alpha_sequence')
            output.append(get_sequence(structure, args.alpha_chain_id, imgt_numbering))

    if args.beta_chain_id:
        for cdr_number, imgt_numbering in enumerate((IMGT_CDR1, IMGT_CDR2, IMGT_CDR3), 1):
            header.append(f'CDR{cdr_number}beta_sequence')
            output.append(get_sequence(structure, args.beta_chain_id, imgt_numbering))

    if args.antigen_chain_id:
        header.append('peptide_sequence')
        output.append(get_sequence(structure, args.antigen_chain_id))

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(','.join(header))
            fh.write('\n')
            fh.write(','.join(output))
            fh.write('\n')

    else:
        print('\t'.join(header))
        print('\t'.join(output))


if __name__ == '__main__':
    main()
