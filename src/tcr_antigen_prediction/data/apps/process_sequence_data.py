"""Process sequence data to create training data for machine learning models.

Required columns:

    - cdr1_alpha
    - cdr2_alpha
    - cdr3_alpha
    - cdr1_beta
    - cdr2_beta
    - cdr3_beta
    - peptide_sequence
    - mhc_pseudo_sequence

"""

import argparse
import logging
import sys

import h5py
import numpy as np
import pandas as pd
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.utils import centre_pad, create_even_folds

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input_data', help='path to the input csv file')
parser.add_argument('--output', '-o', required=True, help='path to output HDF5 file with processed sequences')
parser.add_argument('--seed', default=None, type=int, help='random seed for data splitting')

data_group = parser.add_argument_group('Data')
data_group.add_argument(
    '--negative-proportion',
    default=5,
    type=int,
    help=(
        'ratio of positives to negatives by randomly sampling the background for the negatives (Default: 5, meaning 5 '
        'times more negatives than positives)'
    ),
)

add_logging_arguments(parser)


PROTEIN_LETTERS = sorted(IUPACData.protein_letters_3to1.values())


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    if args.seed:
        logger.info('Setting seed to %d', args.seed)
        rng = np.random.default_rng(args.seed)

    else:
        rng = np.random.default_rng()

    logger.info('Loading data')
    sequence_data = pd.read_csv(args.input_data)
    sequence_data['label'] = [(0, 1)] * len(sequence_data)

    logger.info('Centre padding sequences')
    sequence_data['cdr1_alpha_processed'] = sequence_data['cdr1_alpha'].apply(list).apply(centre_pad, pad_length=8)
    sequence_data['cdr2_alpha_processed'] = sequence_data['cdr2_alpha'].apply(list).apply(centre_pad, pad_length=8)
    sequence_data['cdr3_alpha_processed'] = sequence_data['cdr3_alpha'].apply(list).apply(centre_pad, pad_length=24)

    sequence_data['cdr1_beta_processed'] = sequence_data['cdr1_beta'].apply(list).apply(centre_pad, pad_length=8)
    sequence_data['cdr2_beta_processed'] = sequence_data['cdr2_beta'].apply(list).apply(centre_pad, pad_length=8)
    sequence_data['cdr3_beta_processed'] = sequence_data['cdr3_beta'].apply(list).apply(centre_pad, pad_length=24)

    sequence_data['peptide_processed'] = sequence_data['peptide_sequence'].apply(list).apply(centre_pad, pad_length=12)
    sequence_data['mhc_processed'] = sequence_data['mhc_pseudo_sequence'].apply(list)

    logger.info('One-hot encoding sequences')
    one_hot_mapping = {}

    zero_arr = np.zeros(len(PROTEIN_LETTERS), dtype=int)

    for idx, olc in enumerate(PROTEIN_LETTERS):
        encoding = zero_arr.copy()
        encoding[idx] = 1
        one_hot_mapping[olc] = encoding

    one_hot_mapping['-'] = zero_arr.copy()

    processed_data = sequence_data.filter(regex='_processed$')

    processed_data = processed_data.map(
        lambda seq: np.array([one_hot_mapping[olc] for olc in seq]),
    )

    sequence_data[processed_data.columns] = processed_data

    logger.info('Generating negative data by random sampling')
    sequence_data['collated_cdr_sequences'] = sequence_data.filter(
        regex=r'^cdr[1-3]_(alpha|beta)$',
    ).apply('-'.join, axis='columns')

    peptides = sequence_data['peptide_sequence'].unique()

    tcr_columns = sequence_data.filter(regex='cdr').columns
    pmhc_columns = sequence_data.filter(regex='peptide|mhc').columns

    negative_data = []

    for peptide in peptides:
        in_group = sequence_data[sequence_data['peptide_sequence'] == peptide]
        out_group = sequence_data[
            (sequence_data['peptide_sequence'] != peptide)
            & (~sequence_data['collated_cdr_sequences'].isin(in_group['collated_cdr_sequences']))
        ]

        negatives = (
            out_group.sample(args.negative_proportion * len(in_group), random_state=rng)
            if len(out_group) >= (args.negative_proportion * len(in_group))
            else out_group
        )[tcr_columns]

        negatives = pd.concat(
            [
                negatives,
                pd.DataFrame(
                    [in_group[pmhc_columns].iloc[0].to_numpy()] * len(negatives),
                    columns=in_group[pmhc_columns].iloc[0].index,
                    index=negatives.index,
                ),
            ],
            axis=1,
        )

        negatives['label'] = [(1, 0)] * len(negatives)

        negative_data.append(negatives)

    sequence_data = pd.concat([sequence_data, *negative_data])

    peptide_counts = sequence_data['peptide_sequence'].value_counts()
    peptide_counts = peptide_counts.sort_index().sort_values(ascending=False)
    folds = create_even_folds(
        list(zip(peptide_counts.index.tolist(), peptide_counts.tolist(), strict=True)),
        seed=args.seed,
    )

    sequence_data['fold'] = sequence_data['peptide_sequence'].map(
        {peptide_sequence: i for i, fold in enumerate(folds, 1) for peptide_sequence in fold}
    )

    processed_data = sequence_data.filter(regex='_processed$|label|fold')
    processed_data.columns = [column_name.replace('_processed', '') for column_name in processed_data.columns]

    # TODO standardise column names between apps
    processed_data = processed_data.rename(
        {
            'cdr1_alpha': 'cdr_1a',
            'cdr2_alpha': 'cdr_2a',
            'cdr3_alpha': 'cdr_3a',
            'cdr1_beta': 'cdr_1b',
            'cdr2_beta': 'cdr_2b',
            'cdr3_beta': 'cdr_3b',
            'mhc': 'mhc_pseudo_sequence',
        },
        axis='columns',
    )

    logger.info('Outputting data to %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        for col in processed_data.columns:
            fh[col] = np.array(processed_data[col].tolist())


if __name__ == '__main__':
    main()
