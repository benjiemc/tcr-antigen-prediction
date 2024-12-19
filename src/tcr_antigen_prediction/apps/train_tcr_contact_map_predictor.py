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
from tcr_antigen_prediction.models import TCRContactMapPredictor

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument('training_data', help='path to the training data (HDF5 file)')
inputs.add_argument('--contact-maps', help='path to contact maps for each CDR-peptide/mhc interaction (HDF5 file)')

outputs = parser.add_argument_group('Outputs')
outputs.add_argument('--output', '-o', required=True, help='path to output directory for models')

parser.add_argument('--seed', default=None, type=int, help='Seed for random processes (Default: None)')

data_parameters = parser.add_argument_group('Data Parameters')
data_parameters.add_argument('--cdr-1-length', type=int, default=8, help='maximum CDR 1 length (Default: 8)')
data_parameters.add_argument('--cdr-2-length', type=int, default=8, help='maximum CDR 2 length (Default: 8)')
data_parameters.add_argument('--cdr-3-length', type=int, default=24, help='maximum CDR 3 length (Default: 24)')
data_parameters.add_argument('--peptide-length', type=int, default=12, help='maximum peptide length (Default: 12)')
data_parameters.add_argument(
    '--mhc-pseudo-sequence-length',
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
    cdr_1as: np.ndarray,
    cdr_2as: np.ndarray,
    cdr_3as: np.ndarray,
    cdr_1bs: np.ndarray,
    cdr_2bs: np.ndarray,
    cdr_3bs: np.ndarray,
    peptides: np.ndarray,
    mhc_pseudo_sequences: np.ndarray,
    labels: np.ndarray,
    device: torch.device,
    batch_size: int = 1000,
) -> float:
    """Evaluate model using ROC-AUC.

    Sets model to evaluation mode.
    """
    model.eval()
    indices = np.arange(len(cdr_1as))

    num_eval_batches = int(np.ceil(len(indices) / batch_size))

    predictions = []
    for i in range(num_eval_batches):
        batch_idx_start = i * batch_size

        batch_indices = (
            indices[batch_idx_start : batch_idx_start + batch_size]
            if i < num_eval_batches - 1
            else indices[batch_idx_start:]
        )

        batch_cdr1a = torch.tensor(cdr_1as[batch_indices], dtype=torch.float32, device=device)
        batch_cdr2a = torch.tensor(cdr_2as[batch_indices], dtype=torch.float32, device=device)
        batch_cdr3a = torch.tensor(cdr_3as[batch_indices], dtype=torch.float32, device=device)
        batch_cdr1b = torch.tensor(cdr_1bs[batch_indices], dtype=torch.float32, device=device)
        batch_cdr2b = torch.tensor(cdr_2bs[batch_indices], dtype=torch.float32, device=device)
        batch_cdr3b = torch.tensor(cdr_3bs[batch_indices], dtype=torch.float32, device=device)

        batch_peptide = torch.tensor(peptides[batch_indices], dtype=torch.float32, device=device)
        batch_mhc = torch.tensor(mhc_pseudo_sequences[batch_indices], dtype=torch.float32, device=device)

        prediction = model(
            batch_cdr1a,
            batch_cdr2a,
            batch_cdr3a,
            batch_cdr1b,
            batch_cdr2b,
            batch_cdr3b,
            batch_peptide,
            batch_mhc,
        )

        predictions.append(prediction.detach().cpu().numpy().squeeze(-1))

    predictions = np.concatenate(predictions)

    return roc_auc_score(labels, predictions)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    if args.seed:
        logger.info('Setting random seed to %d', args.seed)
        torch.manual_seed(args.seed)
        random_generator = np.random.default_rng(args.seed)

    else:
        random_generator = np.random.default_rng()

    logger.info('Loading contact maps from %s', args.contact_maps)
    with h5py.File(args.contact_maps) as fh:
        cdr_peptide_contact_maps = (
            torch.tensor(fh['cdr_peptide']['cdr_1a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_peptide']['cdr_2a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_peptide']['cdr_3a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_peptide']['cdr_1b'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_peptide']['cdr_2b'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_peptide']['cdr_3b'][:], dtype=torch.float32),
        )

        cdr_mhc_contact_maps = (
            torch.tensor(fh['cdr_mhc']['cdr_1a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_mhc']['cdr_2a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_mhc']['cdr_3a'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_mhc']['cdr_1b'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_mhc']['cdr_2b'][:], dtype=torch.float32),
            torch.tensor(fh['cdr_mhc']['cdr_3b'][:], dtype=torch.float32),
        )

    logger.info('Loading training data from %s', args.training_data)
    with h5py.File(args.training_data) as fh:
        cdr_1as = fh['cdr_1a'][:]
        cdr_2as = fh['cdr_2a'][:]
        cdr_3as = fh['cdr_3a'][:]
        cdr_1bs = fh['cdr_1b'][:]
        cdr_2bs = fh['cdr_2b'][:]
        cdr_3bs = fh['cdr_3b'][:]

        peptides = fh['peptide'][:]
        mhc_pseudo_sequences = fh['mhc_pseudo_sequence'][:]

        labels = fh['label'][:]
        folds = fh['fold'][:]

    indices = np.arange(len(cdr_1as))

    if not os.path.exists(args.output):
        logger.info('Creating %s', args.output)
        os.mkdir(args.output)

    for fold_idx in np.unique(folds):
        logger.info('Starting fold %d', fold_idx)

        training_indices = indices[folds != fold_idx]
        validation_indices = indices[folds == fold_idx]

        logger.debug('Number of training data points: %d', len(training_indices))
        logger.debug('Number of validation data points: %d', len(validation_indices))

        logger.info('Initialising model')
        model = TCRContactMapPredictor(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            cdr_1_length=args.cdr_1_length,
            cdr_2_length=args.cdr_2_length,
            cdr_3_length=args.cdr_3_length,
            peptide_length=args.peptide_length,
            mhc_length=args.mhc_pseudo_sequence_length,
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

                batch_cdr1a = torch.tensor(cdr_1as[batch_indices], dtype=torch.float32, device=device)
                batch_cdr2a = torch.tensor(cdr_2as[batch_indices], dtype=torch.float32, device=device)
                batch_cdr3a = torch.tensor(cdr_3as[batch_indices], dtype=torch.float32, device=device)
                batch_cdr1b = torch.tensor(cdr_1bs[batch_indices], dtype=torch.float32, device=device)
                batch_cdr2b = torch.tensor(cdr_2bs[batch_indices], dtype=torch.float32, device=device)
                batch_cdr3b = torch.tensor(cdr_3bs[batch_indices], dtype=torch.float32, device=device)

                batch_peptide = torch.tensor(peptides[batch_indices], dtype=torch.float32, device=device)
                batch_mhc = torch.tensor(mhc_pseudo_sequences[batch_indices], dtype=torch.float32, device=device)

                batch_label = torch.tensor(labels[batch_indices], dtype=torch.float32, device=device).unsqueeze(-1)

                prediction = model(
                    batch_cdr1a,
                    batch_cdr2a,
                    batch_cdr3a,
                    batch_cdr1b,
                    batch_cdr2b,
                    batch_cdr3b,
                    batch_peptide,
                    batch_mhc,
                )

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
                        cdr_1as[validation_indices],
                        cdr_2as[validation_indices],
                        cdr_3as[validation_indices],
                        cdr_1bs[validation_indices],
                        cdr_2bs[validation_indices],
                        cdr_3bs[validation_indices],
                        peptides[validation_indices],
                        mhc_pseudo_sequences[validation_indices],
                        labels[validation_indices],
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
