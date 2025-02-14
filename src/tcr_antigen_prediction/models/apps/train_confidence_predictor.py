"""Train a model to provide a confidence score for TCR:pMHC pairing predictions."""

import argparse
import logging
import os
import sys

import h5py
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skl2onnx import to_onnx
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('training_data', help='path to HDF5 file with the training data.')
parser.add_argument('--output', '-o', required=True, help='path to output model file to (ONNX format)')
parser.add_argument(
    '--test-size',
    default=0.15,
    type=float,
    help='Proportion of data to use for evaluation (Default: 0.15)',
)
parser.add_argument('--seed', default=None, type=int, help='Seed for random processes (Default: None)')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    logger.info('Loading data...')
    with h5py.File(args.training_data, 'r') as fh:
        cdr_1as = fh['cdr1_alpha'][:]
        cdr_2as = fh['cdr2_alpha'][:]
        cdr_3as = fh['cdr3_alpha'][:]
        cdr_1bs = fh['cdr1_beta'][:]
        cdr_2bs = fh['cdr2_beta'][:]
        cdr_3bs = fh['cdr3_beta'][:]
        peptides = fh['peptide'][:]
        mhcs = fh['mhc_pseudo'][:]

        folds = fh['fold'][:]

    indices = np.arange(len(cdr_1as))

    tcr_pmhc = np.concatenate([cdr_1as, cdr_2as, cdr_3as, cdr_1bs, cdr_2bs, cdr_3bs, peptides, mhcs], axis=1)
    tcr_pmhc_flat = tcr_pmhc.reshape(tcr_pmhc.shape[0], -1)

    logger.info('Computing distances')
    tcr_pmhc_distances = squareform(pdist(tcr_pmhc_flat, metric='euclidean'))

    logger.debug('Splitting training and testing points')
    train_indices, test_indices = train_test_split(indices, test_size=args.test_size, random_state=args.seed)

    logger.debug('Train indices: %s', list(train_indices))
    logger.info('Test indices: %s', list(test_indices))

    for fold in np.unique(folds):
        logger.info('Starting model for fold %d', fold)
        training_idxs = indices[folds != fold]

        logger.debug('Calculating average distance to training data')
        tcr_pmhc_distances_avg = np.mean(tcr_pmhc_distances[:, training_idxs], axis=1)

        logger.debug('Normalising distances')
        tcr_pmhc_distances_avg_norm = -tcr_pmhc_distances_avg

        min_val = np.min(tcr_pmhc_distances_avg_norm)
        max_val = np.max(tcr_pmhc_distances_avg_norm)

        tcr_pmhc_distances_avg_norm = (tcr_pmhc_distances_avg_norm - min_val) / (max_val - min_val)

        logger.info('Training confidence score model')
        confidence_model = RandomForestRegressor(random_state=args.seed).fit(
            tcr_pmhc_flat[train_indices],
            tcr_pmhc_distances_avg_norm[train_indices],
        )

        logger.info('Evaluating model')
        r_squared = confidence_model.score(tcr_pmhc_flat[test_indices], tcr_pmhc_distances_avg_norm[test_indices])
        logger.info('R^2: %f', r_squared)

        predictions = confidence_model.predict(tcr_pmhc_flat[test_indices])
        rmse = np.sqrt(np.mean((tcr_pmhc_distances_avg_norm[test_indices] - predictions) ** 2))
        logger.info('RMSE: %f', rmse)

        logger.info('Saving model')
        onnx = to_onnx(confidence_model, tcr_pmhc_flat[:1].astype(np.float32))

        output_name = os.path.join(
            os.path.dirname(args.output),
            f'_{fold}.'.join(os.path.basename(args.output).split('.', 1)),
        )

        with open(output_name, 'wb') as fh:
            fh.write(onnx.SerializeToString())


if __name__ == '__main__':
    main()
