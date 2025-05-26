"""Process sequence data to create training data for machine learning models.

Required columns:

    - cdr1_alpha
    - cdr2_alpha
    - cdr3_alpha
    - cdr1_beta
    - cdr2_beta
    - cdr3_beta
    - peptide
    - mhc_pseudo
    - label
    - fold

"""

import argparse
import logging
import sys

import h5py
import numpy as np
import pandas as pd

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.amino_acid_encodings import BLOSUM_50_ENCODING, ONE_HOT_ENCODING
from tcr_antigen_prediction.data.utils import centre_pad, left_pad, right_pad

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input_data', help='path to the input csv file')
parser.add_argument('--output', '-o', required=True, help='path to output HDF5 file with processed sequences')

data_group = parser.add_argument_group('Data')
data_group.add_argument(
    '--cdr1-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 1a sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr2-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 2a sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr3-alpha-length',
    default=24,
    type=int,
    help='Length to pad CDR 3a sequences to (Default: 24)',
)
data_group.add_argument(
    '--cdr1-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 1b sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr2-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 2b sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr3-beta-length',
    default=24,
    type=int,
    help='Length to pad CDR 3b sequences to (Default: 24)',
)
data_group.add_argument(
    '--peptide-length',
    default=12,
    type=int,
    help='Length to pad peptide sequences to (Default: 12)',
)
data_group.add_argument(
    '--pad-direction',
    choices=['centre', 'left', 'right'],
    default='centre',
    help="Direction to pad sequences (Default: 'centre')",
)
data_group.add_argument(
    '--encoding',
    choices=['blosum50', 'one-hot'],
    default='one-hot',
    help='numerical encoding to use for the sequences.',
)
data_group.add_argument('--normalisation-factor', type=float, default=1.0, help='normalise encoding values by factor')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading data')
    sequence_data = pd.read_csv(args.input_data)

    match args.pad_direction:
        case 'centre':
            logger.info('Centre padding sequences')
            pad_func = centre_pad

        case 'right':
            logger.info('Right padding sequences')
            pad_func = right_pad

        case 'left':
            logger.info('Left padding sequences')
            pad_func = left_pad

    sequence_data['cdr1_alpha_processed'] = (
        sequence_data['cdr1_alpha'].apply(list).apply(pad_func, pad_length=args.cdr1_alpha_length)
    )
    sequence_data['cdr2_alpha_processed'] = (
        sequence_data['cdr2_alpha'].apply(list).apply(pad_func, pad_length=args.cdr2_alpha_length)
    )
    sequence_data['cdr3_alpha_processed'] = (
        sequence_data['cdr3_alpha'].apply(list).apply(pad_func, pad_length=args.cdr3_alpha_length)
    )
    sequence_data['cdr1_beta_processed'] = (
        sequence_data['cdr1_beta'].apply(list).apply(pad_func, pad_length=args.cdr1_beta_length)
    )
    sequence_data['cdr2_beta_processed'] = (
        sequence_data['cdr2_beta'].apply(list).apply(pad_func, pad_length=args.cdr2_beta_length)
    )
    sequence_data['cdr3_beta_processed'] = (
        sequence_data['cdr3_beta'].apply(list).apply(pad_func, pad_length=args.cdr3_beta_length)
    )

    sequence_data['peptide_processed'] = (
        sequence_data['peptide'].apply(list).apply(pad_func, pad_length=args.peptide_length)
    )
    sequence_data['mhc_pseudo_processed'] = sequence_data['mhc_pseudo'].apply(list)

    match args.encoding:
        case 'one-hot':
            logger.info('One-hot encoding sequences')
            encoding = ONE_HOT_ENCODING

        case 'blosum50':
            logger.info('Blosum 50 encoding sequences')
            encoding = BLOSUM_50_ENCODING

    processed_data = sequence_data.filter(regex='_processed$')
    processed_data = processed_data.map(
        lambda seq: np.array([encoding[olc] for olc in seq]),
    )
    logger.debug('Normalising by factor %f', args.normalisation_factor)
    processed_data = processed_data.div(args.normalisation_factor)

    sequence_data[processed_data.columns] = processed_data

    processed_data = sequence_data.filter(regex='_processed$|label|fold')
    processed_data.columns = [column_name.replace('_processed', '') for column_name in processed_data.columns]

    logger.info('Outputting data to %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        for col in processed_data.columns:
            fh[col] = np.array(processed_data[col].tolist())


if __name__ == '__main__':
    main()
