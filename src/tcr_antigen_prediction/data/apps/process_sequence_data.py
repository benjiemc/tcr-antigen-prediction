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

"""

import argparse
import logging
import sys

import h5py
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.amino_acid_encodings import BLOSUM_50_ENCODING, ONE_HOT_ENCODING
from tcr_antigen_prediction.data.utils import centre_pad, create_even_folds, find_common_groups, left_pad, right_pad

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
data_group.add_argument(
    '--num-folds',
    type=int,
    default=5,
    help='Number of folds to create in the dataset (Default: 5)',
)
data_group.add_argument(
    '--split-type',
    choices=['random', 'tcr', 'peptide', 'pMHC', 'levenshtein'],
    default='peptide',
    help=(
        "Method to partition data between folds (Default: 'peptide'). 'random' means to randomly shuffle data between"
        " folds, 'tcr' means no TCRs are shared across folds, 'peptide' means to ensure no peptides are shared across "
        "folds, 'pMHC' means that no peptides or MHCs (based on pseudo sequence) are shared across folds, and "
        "'levenshtein' means to use a levenshtein distance to separate data points between folds (more parameters "
        "below)."
    ),
)
data_group.add_argument(
    '--split-distance',
    default=6,
    type=int,
    help=(
        "If using a 'levenshtein' split, the mimimum distance between data points in two different folds (Default: 6)."
    ),
)
data_group.add_argument(
    '--split-entities',
    nargs='+',
    choices=[
        'cdr1_alpha',
        'cdr2_alpha',
        'cdr3_alpha',
        'cdr1_beta',
        'cdr2_beta',
        'cdr3_beta',
        'peptide',
        'mhc_pseudo',
    ],
    help="If using a 'levenshtein' split, the entities to consider as part of the split.",
)
data_group.add_argument(
    '--distances',
    help=(
        "If using a 'levenshtein' split, the path to an hdf5 file containing the distances between each type of entity "
        "considered."
    ),
)

add_logging_arguments(parser)


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
    sequence_data['label'] = 1

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

    logger.info('Generating negative data by random sampling')
    sequence_data['collated_cdrs'] = sequence_data.filter(
        regex=r'^cdr[1-3]_(alpha|beta)$',
    ).apply('-'.join, axis='columns')

    peptides = sequence_data['peptide'].unique()

    tcr_columns = sequence_data.filter(regex='cdr').columns
    pmhc_columns = sequence_data.filter(regex='peptide|mhc').columns

    negative_data = []

    for peptide in peptides:
        in_group = sequence_data[sequence_data['peptide'] == peptide]
        out_group = sequence_data[
            (sequence_data['peptide'] != peptide) & (~sequence_data['collated_cdrs'].isin(in_group['collated_cdrs']))
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

        negatives['label'] = 0

        negative_data.append(negatives)

    sequence_data = pd.concat([sequence_data, *negative_data]).reset_index(drop=True)

    logger.info('Creating %d cross-validation folds', args.num_folds)
    match args.split_type:
        case 'random':
            logger.debug('Randomly assigning cross-validation folds')

            if len(sequence_data) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d folds will be created',
                    args.num_folds,
                    len(sequence_data),
                )

            folds = create_even_folds(
                [(idx, 1) for idx in sequence_data.index],
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data.index.map({idx: i for i, fold in enumerate(folds, 1) for idx in fold})

        case 'tcr':
            logger.debug('Splitting TCRs across cross-validation folds')
            tcr_counts = sequence_data['collated_cdrs'].value_counts()
            tcr_counts = tcr_counts.sort_index().sort_values(ascending=False)

            if len(tcr_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(tcr_counts),
                )

            folds = create_even_folds(
                list(zip(tcr_counts.index.tolist(), tcr_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data['collated_cdrs'].map(
                {tcr_sequence: i for i, fold in enumerate(folds, 1) for tcr_sequence in fold}
            )

        case 'peptide':
            logger.debug('Splitting peptides across cross-validation folds')
            peptide_counts = sequence_data['peptide'].value_counts()
            peptide_counts = peptide_counts.sort_index().sort_values(ascending=False)

            if len(peptide_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(peptide_counts),
                )

            folds = create_even_folds(
                list(zip(peptide_counts.index.tolist(), peptide_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data['peptide'].map(
                {peptide_sequence: i for i, fold in enumerate(folds, 1) for peptide_sequence in fold}
            )

        case 'pMHC':
            logger.debug('Splitting peptides across cross-validation folds')
            pmhc_groups = pd.Series(find_common_groups(sequence_data, ['peptide', 'mhc_pseudo']))
            pmhc_group_counts = pmhc_groups.value_counts()

            if len(pmhc_group_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(pmhc_group_counts),
                )

            folds = create_even_folds(
                list(zip(pmhc_group_counts.index.tolist(), pmhc_group_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = pmhc_groups.map(
                {group_id: i for i, fold in enumerate(folds, 1) for group_id in fold}
            )

        case 'levenshtein':
            logger.debug(
                'Splitting data points across cross-validation folds based on a levenshtein distance of %d',
                args.split_distance,
            )

            with h5py.File(args.distances) as fh:
                distances = np.zeros((len(sequence_data), len(sequence_data)), dtype=int)

                for entity in args.split_entities:
                    names_map = {name.decode('utf-8'): idx for idx, name in enumerate(fh[entity]['names'][:])}
                    indices = sequence_data[entity].map(names_map).to_numpy()

                    distances += fh[entity]['distance_matrix'][:][np.ix_(indices, indices)]

            clusters = AgglomerativeClustering(
                n_clusters=None,
                distance_threshold=args.split_distance,
                metric='precomputed',
                linkage='single',
            ).fit_predict(distances)

            cluster_sizes = {}
            for cluster in clusters:
                if cluster in cluster_sizes:
                    cluster_sizes[cluster] += 1

                else:
                    cluster_sizes[cluster] = 1

            cluster_sizes = sorted(cluster_sizes.items(), key=lambda cluster_size: cluster_size[1], reverse=True)

            if len(cluster_sizes) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(cluster_sizes),
                )

            folds = create_even_folds(cluster_sizes, num_folds=args.num_folds, seed=args.seed)
            fold_map = {cluster: i for i, fold in enumerate(folds, 1) for cluster in fold}
            sequence_data['fold'] = [fold_map[cluster] for cluster in clusters]

    processed_data = sequence_data.filter(regex='_processed$|label|fold')
    processed_data.columns = [column_name.replace('_processed', '') for column_name in processed_data.columns]

    logger.info('Outputting data to %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        for col in processed_data.columns:
            fh[col] = np.array(processed_data[col].tolist())


if __name__ == '__main__':
    main()
