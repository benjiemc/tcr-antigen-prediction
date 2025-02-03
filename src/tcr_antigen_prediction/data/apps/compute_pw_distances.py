"""Compute pair-wise levenshtein distances between input strings."""

import argparse
import logging
import sys

import numpy as np
from pyxdameraulevenshtein import damerau_levenshtein_distance

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input', nargs='+', help='input strings')
parser.add_argument('--output', '-o', help='path to output txt file')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.debug('Computing distances')
    distance_matrix = np.zeros((len(args.input), len(args.input)), dtype=np.int32)

    for i, seq1 in enumerate(args.input):
        for j, seq2 in enumerate(args.input[i + 1 :]):
            distance_matrix[i, j + i + 1] = damerau_levenshtein_distance(seq1, seq2)

    logger.debug('Symmetrising distance matrix')
    distance_matrix = np.maximum(distance_matrix, distance_matrix.T)

    logger.debug('Ouputting matrix to %s', args.output)
    np.savetxt(args.output, distance_matrix, fmt='%i')


if __name__ == '__main__':
    main()
