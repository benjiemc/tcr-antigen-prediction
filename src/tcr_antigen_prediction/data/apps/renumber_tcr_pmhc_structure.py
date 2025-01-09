"""Renumber TCR structure following either IMGT or Aho numbering.

Requirements:
    - ANARCI: https://github.com/oxpig/ANARCI

"""

import argparse
import logging
import sys

from Bio.PDB import PDBIO, Model, PDBParser, Structure

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import renumber_chain
from tcr_antigen_prediction.data.structure import get_header

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to the pdb structure file')
parser.add_argument('--output', '-o', help='name of output structure file')

add_logging_arguments(parser)


def main():
    """Entry point for script."""
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    with open(args.structure, 'r') as fh:
        header = get_header(fh.read())
        fh.seek(0)

        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('', fh)

    chain_map = []
    output_structure = Structure.Structure('otuput')
    for model in structure:
        new_model = Model.Model(model.id)

        for chain in model:
            try:
                renumbered_chain, chain_type = renumber_chain(chain)

            except ValueError:
                logger.debug('Chain ID %s not identified by ANARCI', chain.id)
                new_model.add(chain.copy())

            else:
                chain_map.append((chain_type, chain.id))
                new_model.add(renumbered_chain)

        output_structure.add(new_model)

    with open(args.output, 'w') as fh:
        fh.write('REMARK     Renumbered using IMGT numbering provided by ANARCI (DOI: 10.1093/bioinformatics/btv552)')
        fh.write('\n')
        fh.write(f"REMARK     {' '.join([f'{chain_type}CHAIN={chain_id}' for chain_type, chain_id in chain_map])}")
        fh.write('\n')

        if len(header) > 0:
            fh.write(header)
            fh.write('\n')

        io = PDBIO()
        io.set_structure(output_structure)
        io.save(fh)


if __name__ == '__main__':
    main()
