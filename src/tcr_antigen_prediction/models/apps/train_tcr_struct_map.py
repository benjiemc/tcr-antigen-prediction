"""Train TCRContactMapPredictor model."""

import argparse
import logging
import os
import sys

import h5py
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from torch import nn, optim

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.models import TCRStructMap

logger = logging.getLogger()

MODEL_FEATURES = [
    'cdr1_alpha',
    'cdr2_alpha',
    'cdr3_alpha',
    'cdr1_beta',
    'cdr2_beta',
    'cdr3_beta',
    'peptide',
    'mhc_pseudo',
]

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument('training_data', help='path to the training data (HDF5 file)')
inputs.add_argument('--contact-maps', help='path to contact maps for each CDR-peptide/mhc interaction (HDF5 file)')
inputs.add_argument('--contact-probabilities', help='path to table of amino acid pairing frequencies')

outputs = parser.add_argument_group('Outputs')
outputs.add_argument('--output', '-o', required=True, help='path to output directory for models')

parser.add_argument('--seed', default=None, type=int, help='Seed for random processes (Default: None)')

data_parameters = parser.add_argument_group('Data Parameters')
data_parameters.add_argument(
    '--features-to-include',
    nargs='+',
    default='all',
    choices=[*MODEL_FEATURES, 'all'],
    help=(
        'features to include in training (Default: all). This is useful for understanding the impact of different '
        'features on model perfomance.'
    ),
)
data_parameters.add_argument('--cdr1-alpha-length', type=int, default=8, help='maximum CDR1-alpha length (Default: 8)')
data_parameters.add_argument('--cdr2-alpha-length', type=int, default=8, help='maximum CDR2-alpha length (Default: 8)')
data_parameters.add_argument(
    '--cdr3-alpha-length',
    type=int,
    default=24,
    help='maximum CDR3-alpha length (Default: 24)',
)
data_parameters.add_argument('--cdr1-beta-length', type=int, default=8, help='maximum CDR1-beta length (Default: 8)')
data_parameters.add_argument('--cdr2-beta-length', type=int, default=8, help='maximum CDR2-beta length (Default: 8)')
data_parameters.add_argument('--cdr3-beta-length', type=int, default=24, help='maximum CDR3-beta length (Default: 24)')
data_parameters.add_argument('--peptide-length', type=int, default=12, help='maximum peptide length (Default: 12)')
data_parameters.add_argument(
    '--mhc-pseudo-length',
    type=int,
    default=26,
    help='length of MHC pseudo sequences (Default: 26)',
)

training_params = parser.add_argument_group('Training Parameters')
training_params.add_argument('--batch-size', type=int, default=32, help='training batch size (Default: 32)')
training_params.add_argument('--num-epochs', type=int, default=25, help='number of epochs (Default: 25)')
training_params.add_argument(
    '--learning-rate',
    type=float,
    default=1e-3,
    help='optimizer learning rate (Default: 1e-3)',
)
training_params.add_argument(
    '--drop-out-rate',
    default=0.2,
    type=float,
    help='Drop-out rate for training (Default: 0.2)',
)
training_params.add_argument(
    '--learn-contact-maps',
    action='store_true',
    help='Unfreeze contact map values during training',
)
training_params.add_argument(
    '--checkpoint-interval', type=int, default=1000, help='interval to save the model during training (Default: 1000)'
)

validation_params = parser.add_argument_group('Validation Parameters')
validation_params.add_argument(
    '--eval-interval',
    type=int,
    default=500,
    help='Interval to evaluate model performance (Default: 500)',
)
validation_params.add_argument(
    '--eval-batch-size', type=int, default=1000, help='batch size for evaluating the model (Default: 1000)'
)

add_logging_arguments(parser)


def evaluate_model(
    model: nn.Module,
    input_data: dict[str, np.ndarray],
    labels: np.ndarray,
    device: torch.device,
    batch_size: int = 1000,
) -> float:
    """Evaluate model using ROC-AUC.

    Sets model to evaluation mode.
    """
    model.eval()
    indices = np.arange(len(labels))

    num_eval_batches = int(np.ceil(len(indices) / batch_size))

    predictions = []
    for i in range(num_eval_batches):
        batch_idx_start = i * batch_size

        batch_indices = (
            indices[batch_idx_start : batch_idx_start + batch_size]
            if i < num_eval_batches - 1
            else indices[batch_idx_start:]
        )

        batch = {
            feature: torch.tensor(data[batch_indices], dtype=torch.float32, device=device)
            for feature, data in input_data.items()
        }

        prediction = model(**batch)

        predictions.append(prediction.detach().cpu().numpy().squeeze(-1))

    predictions = np.concatenate(predictions)

    return roc_auc_score(labels, predictions)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    feature_order = {name: num for num, name in enumerate(MODEL_FEATURES)}
    features = (
        MODEL_FEATURES
        if args.features_to_include == 'all'
        else sorted(args.features_to_include, key=lambda name: feature_order[name])
    )

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    if args.seed:
        logger.info('Setting random seed to %d', args.seed)
        torch.manual_seed(args.seed)
        random_generator = np.random.default_rng(args.seed)

    else:
        random_generator = np.random.default_rng()

    if args.contact_maps:
        logger.info('Loading contact maps from %s', args.contact_maps)
        with h5py.File(args.contact_maps) as fh:
            cdr_peptide_contact_maps = (
                torch.tensor(fh['peptide']['cdr1_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['peptide']['cdr2_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['peptide']['cdr3_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['peptide']['cdr1_beta'][:], dtype=torch.float32),
                torch.tensor(fh['peptide']['cdr2_beta'][:], dtype=torch.float32),
                torch.tensor(fh['peptide']['cdr3_beta'][:], dtype=torch.float32),
            )

            cdr_mhc_contact_maps = (
                torch.tensor(fh['mhc_pseudo']['cdr1_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['mhc_pseudo']['cdr2_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['mhc_pseudo']['cdr3_alpha'][:], dtype=torch.float32),
                torch.tensor(fh['mhc_pseudo']['cdr1_beta'][:], dtype=torch.float32),
                torch.tensor(fh['mhc_pseudo']['cdr2_beta'][:], dtype=torch.float32),
                torch.tensor(fh['mhc_pseudo']['cdr3_beta'][:], dtype=torch.float32),
            )

    else:
        cdr_peptide_contact_maps = None
        cdr_mhc_contact_maps = None

    if args.contact_probabilities:
        logger.info('Loading contact probabilities from %s', args.contact_probabilities)
        contact_probabilities = torch.tensor(np.loadtxt(args.contact_probabilities), dtype=torch.float32)

    else:
        contact_probabilities = None

    logger.info('Loading training data from %s', args.training_data)
    with h5py.File(args.training_data) as fh:
        input_data = {name: fh[name][:] for name in features}

        labels = fh['label'][:]
        folds = fh['fold'][:]

    indices = np.arange(len(labels))

    if not os.path.exists(args.output):
        logger.info('Creating %s', args.output)
        os.mkdir(args.output)

    for fold_idx in np.unique(folds):
        logger.info('Starting fold %d', fold_idx)

        training_indices = indices[folds != fold_idx]
        validation_indices = indices[folds == fold_idx]

        validation_data = {name: data[validation_indices] for name, data in input_data.items()}
        validation_labels = labels[validation_indices]

        logger.debug('Number of training data points: %d', len(training_indices))
        logger.debug('Number of validation data points: %d', len(validation_indices))

        logger.info('Initialising model')
        model = TCRStructMap(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            contact_probabilities,
            cdr1_alpha_length=args.cdr1_alpha_length if 'cdr1_alpha' in features else 0,
            cdr2_alpha_length=args.cdr2_alpha_length if 'cdr2_alpha' in features else 0,
            cdr3_alpha_length=args.cdr3_alpha_length if 'cdr3_alpha' in features else 0,
            cdr1_beta_length=args.cdr1_beta_length if 'cdr1_beta' in features else 0,
            cdr2_beta_length=args.cdr2_beta_length if 'cdr2_beta' in features else 0,
            cdr3_beta_length=args.cdr3_beta_length if 'cdr3_beta' in features else 0,
            peptide_length=args.peptide_length if 'peptide' in features else 0,
            mhc_length=args.mhc_pseudo_length if 'mhc_pseudo' in features else 0,
            drop_out_rate=args.drop_out_rate,
            learn_contact_maps=args.learn_contact_maps,
        )
        model.train()

        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        logger.debug('Transferring model to %s', str(device))
        model.to(device)

        loss_fn = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=args.learning_rate)

        output_name = os.path.join(args.output, f'model_{fold_idx}.pt')

        batches_per_epoch = len(training_indices) // args.batch_size

        logger.info('Starting training')
        for epoch in range(1, args.num_epochs + 1):
            logger.info('Starting epoch %d', epoch)
            epoch_indices = training_indices.copy()
            random_generator.shuffle(epoch_indices)

            for i in range(batches_per_epoch):
                logger.debug('Training batch %d of %d', i + 1, batches_per_epoch)
                batch_indices = epoch_indices[i * args.batch_size : i * args.batch_size + args.batch_size]

                batch = {
                    feature: torch.tensor(data[batch_indices], dtype=torch.float32, device=device)
                    for feature, data in input_data.items()
                }

                batch_label = torch.tensor(labels[batch_indices], dtype=torch.float32, device=device).unsqueeze(-1)

                prediction = model(**batch)

                loss = loss_fn(prediction, batch_label)
                loss.backward()

                optimizer.step()
                optimizer.zero_grad()

                if i % args.eval_interval == 0:
                    logger.info('Evaluating at iteration %d', (i + 1) + (batches_per_epoch * (epoch - 1)))

                    training_roc_auc = roc_auc_score(
                        batch_label.detach().cpu().numpy().squeeze(-1),
                        prediction.detach().cpu().numpy().squeeze(-1),
                    )

                    validation_roc_auc = evaluate_model(
                        model,
                        validation_data,
                        validation_labels,
                        device,
                        batch_size=args.eval_batch_size,
                    )

                    logger.info('Loss: %f', loss.item())
                    logger.info('Training ROC-AUC: %f', training_roc_auc)
                    logger.info('Validation ROC-AUC: %f', validation_roc_auc)

                    model.train()

                if i % args.checkpoint_interval == 0:
                    logger.debug('CHECKPOINT: Saving model to %s', output_name)
                    torch.save(model.state_dict(), output_name)

        logger.info('Saving model to %s', output_name)
        torch.save(model.state_dict(), output_name)


if __name__ == '__main__':
    main()
