"""Replace HETATM records with ATOM records if they are valid residues."""

import argparse
import logging
import sys

from Bio.PDB import PDBIO, PDBParser

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.structure import PROTEIN_LETTERS, get_header

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to pdb file')
parser.add_argument('--output', '-o', help='path to output pdb file')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    with open(args.structure, 'r') as fh:
        header = get_header(fh.read())
        fh.seek(0)

        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('', fh)

    for model in structure:
        for chain in model:
            for res in chain:
                # Does not account for non-standard residues
                if res.id[0] == ' ' or res.get_resname() not in PROTEIN_LETTERS:
                    continue

                res.id = (' ', res.id[1], res.id[2])

    io = PDBIO()
    io.set_structure(structure)

    fh = open(args.output, 'w') if args.output else sys.stdout  # noqa: SIM115

    fh.write(header)
    io.save(fh)

    if args.output:
        fh.close()


if __name__ == '__main__':
    main()
