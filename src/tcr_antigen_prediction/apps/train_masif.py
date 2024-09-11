# This file is part of the MaSIF project.
#
# Copyright 2019 - Gainza P, Sverrisson F, Monti F, Rodola, Bronstein MM, Correia BE
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
'''Train or finetune MaSIF-PPI model.

TODO Improve generalisability of module.

'''
import argparse
import logging
import os
from typing import List, Optional, Tuple

import numpy as np
import pymesh
from sklearn import metrics
from scipy.spatial import cKDTree

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.models import MasifPPISearch


logger = logging.getLogger()

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('training', nargs='+', help='path to training data')
parser.add_argument('--output', '-o', required=True, help='path to model output')
parser.add_argument('--validation-data', nargs='+', help='path to validation data')
parser.add_argument('--ply-dir', help='path to ply mesh files')
parser.add_argument('--model', help='path to the model for restarting training or fine tuning.')
parser.add_argument('--seed', default=None, type=int, help='seed for random state')
parser.add_argument('--binder-name', required=True, help='name of binders')
parser.add_argument('--positive-name', required=True, help='name of positives')

training_parameters = parser.add_argument_group('Training Parameters')
training_parameters.add_argument('--num-iterations', type=int, default=1_000_000,
                                 help='number of training iterations (Default: 1_000_000)')
training_parameters.add_argument('--num-iter-eval', type=int, default=1000,
                                 help='number of tranining loops before evaluating the model (Default: 1000)')
training_parameters.add_argument('--batch-size', type=int, default=32, help='batch size (Default: 32)')
training_parameters.add_argument('--validation-batch-size', type=int, default=1000,
                                 help='batch size of validation/testing (Default: 1000)')

data_parameters = parser.add_argument_group('Data Parameters')
data_parameters.add_argument('--feat-mask', nargs='+', type=float, default=[1.0, 1.0, 1.0, 1.0, 1.0],
                             help='feature mask (Default: [1.0, 1.0, 1.0, 1.0, 1.0])')
data_parameters.add_argument('--max-distance', type=float, default=12.0,
                             help='radius for the neural network (Default: 12.0)')
data_parameters.add_argument('--contact-distance', type=float, default=5.0,
                             help='cutoff to consider patches in contact (Default: 5.0 Å)')
data_parameters.add_argument('--sc-max-cutoff', type=float, default=None,
                             help='limit on shape complementarity upper bound (Default: None)')
data_parameters.add_argument('--sc-min-cutoff', type=float, default=None,
                             help='limit on shape complementarity lower bound (Default: None)')
data_parameters.add_argument('--pos-surf-accept-probability', type=float, default=1.0,
                             help='TODO (Default: 1.0)')

add_logging_arguments(parser)


def mask_input_feat(input_feat, mask):
    '''Apply mask to input_feat'''
    mymask = np.where(np.array(mask) == 0.0)[0]
    return np.delete(input_feat, mymask, axis=2)


def construct_batch(binder_rho_wrt_center,
                    binder_theta_wrt_center,
                    binder_input_feat,
                    binder_mask,
                    c_pos_training_idx,
                    pos_rho_wrt_center,
                    pos_theta_wrt_center,
                    pos_input_feat,
                    pos_mask,
                    c_neg_training_idx,
                    neg_rho_wrt_center,
                    neg_theta_wrt_center,
                    neg_input_feat,
                    neg_mask):
    '''
    Construct a batch of training data. Features and theta are flipped for the binder in construct_batch (except for
    hydrophobicity).
    '''
    batch_rho_coords_binder = np.expand_dims(
        binder_rho_wrt_center[c_pos_training_idx], 2
    )
    batch_theta_coords_binder = np.expand_dims(
        binder_theta_wrt_center[c_pos_training_idx], 2
    )
    batch_input_feat_binder = binder_input_feat[c_pos_training_idx]
    batch_mask_binder = binder_mask[c_pos_training_idx]

    batch_rho_coords_pos = np.expand_dims(pos_rho_wrt_center[c_pos_training_idx], 2)
    batch_theta_coords_pos = np.expand_dims(pos_theta_wrt_center[c_pos_training_idx], 2)
    batch_input_feat_pos = pos_input_feat[c_pos_training_idx]
    batch_mask_pos = pos_mask[c_pos_training_idx]

    # Negate the input_features of the binder, except the last column.
    batch_input_feat_binder = -batch_input_feat_binder
    # TODO: This should not be like this ... it is a hack.
    if batch_input_feat_binder.shape[2] == 5 or batch_input_feat_binder.shape[2] == 3:
        batch_input_feat_binder[:, :, -1] = -batch_input_feat_binder[
            :, :, -1
        ]  # Do not negate hydrophobicity.
    # Also negate the theta coords for the binder.
    batch_theta_coords_binder = 2 * np.pi - batch_theta_coords_binder

    batch_rho_coords_neg = np.expand_dims(neg_rho_wrt_center[c_neg_training_idx], 2)
    batch_theta_coords_neg = np.expand_dims(neg_theta_wrt_center[c_neg_training_idx], 2)
    batch_input_feat_neg = neg_input_feat[c_neg_training_idx]
    batch_mask_neg = neg_mask[c_neg_training_idx]

    batch_rho_coords_neg_2 = batch_rho_coords_binder.copy()
    batch_theta_coords_neg_2 = batch_theta_coords_binder.copy()
    batch_input_feat_neg_2 = batch_input_feat_binder.copy()
    batch_mask_neg_2 = batch_mask_binder.copy()

    batch_rho_coords = np.concatenate(
        [
            batch_rho_coords_pos,
            batch_rho_coords_binder,
            batch_rho_coords_neg,
            batch_rho_coords_neg_2,
        ],
        axis=0,
    )
    batch_theta_coords = np.concatenate(
        [
            batch_theta_coords_pos,
            batch_theta_coords_binder,
            batch_theta_coords_neg,
            batch_theta_coords_neg_2,
        ],
        axis=0,
    )
    batch_input_feat = np.concatenate(
        [
            batch_input_feat_pos,
            batch_input_feat_binder,
            batch_input_feat_neg,
            batch_input_feat_neg_2,
        ],
        axis=0,
    )
    batch_mask = np.concatenate(
        [batch_mask_pos, batch_mask_binder, batch_mask_neg, batch_mask_neg_2], axis=0
    )
    # expand the last dimension of the mask (batch_size, max_points_patch, 1)
    batch_mask = np.expand_dims(batch_mask, 2)

    return batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask


def compute_dists(descs1, descs2):
    dists = np.sqrt(np.sum(np.square(descs1 - descs2), axis=1))
    return dists


def construct_batch_val_test(c_idx, rho_wrt_center, theta_wrt_center, input_feat, mask, flip=False):
    batch_rho_coords = np.expand_dims(rho_wrt_center[c_idx], 2)
    batch_theta_coords = np.expand_dims(theta_wrt_center[c_idx], 2)
    batch_input_feat = input_feat[c_idx]
    batch_mask = mask[c_idx]
    batch_mask = np.expand_dims(batch_mask, 2)
    # Flip features and theta (except hydrophobicity)
    if flip:
        batch_input_feat = -batch_input_feat
        batch_theta_coords = 2 * np.pi - batch_theta_coords
        assert len(batch_input_feat.shape) == 3
        # Hydrophobicity is not flipped. -- Fix this.
        if batch_input_feat.shape[2] == 5 or batch_input_feat.shape[2] == 3:
            batch_input_feat[:, :, -1] = -batch_input_feat[:, :, -1]

    return batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask


def compute_val_test_desc(learning_obj,
                          idx,
                          rho_wrt_center,
                          theta_wrt_center,
                          input_feat,
                          mask,
                          batch_size=100,
                          flip=False):
    all_descs = []
    num_batches = int(np.ceil(float(len(idx)) / float(batch_size)))

    # Compute all desc for positive shapes.
    for kk in range(num_batches):
        c_idx = idx[np.arange(kk * batch_size, min((kk + 1) * batch_size, len(idx)))]

        batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask = construct_batch_val_test(c_idx,
                                                                                                      rho_wrt_center,
                                                                                                      theta_wrt_center,
                                                                                                      input_feat,
                                                                                                      mask,
                                                                                                      flip=flip)

        feed_dict = {
            learning_obj.rho_coords: batch_rho_coords,
            learning_obj.theta_coords: batch_theta_coords,
            learning_obj.input_feat: batch_input_feat,
            learning_obj.mask: batch_mask,
            learning_obj.keep_prob: 1.0,
        }

        desc = learning_obj.session.run([learning_obj.global_desc], feed_dict=feed_dict)
        desc = np.squeeze(desc)

        if len(desc.shape) == 1:
            desc = np.expand_dims(desc, 0)

        all_descs.append(desc)

    if len(all_descs) > 1:
        all_descs = np.concatenate(all_descs, axis=0)

    else:
        all_descs = all_descs[0]

    return all_descs


def compute_roc_auc(pos, neg):
    labels = np.concatenate([np.ones((len(pos))), np.zeros((len(neg)))])
    dist_pairs = np.concatenate([pos, neg])
    return metrics.roc_auc_score(labels, dist_pairs)


def aggregate_data(data_dirs: List[str],
                   binder_name: str,
                   positive_name: str,
                   ply_dir: str,
                   contact_distance: float = 5.0,
                   pos_surf_accept_probability: float = 1.0,
                   sc_min_filt: Optional[float] = None,
                   sc_max_filt: Optional[float] = None) -> Tuple:
    binder_rho_wrt_center = []
    binder_theta_wrt_center = []
    binder_input_feat = []
    binder_mask = []

    pos_rho_wrt_center = []
    pos_theta_wrt_center = []
    pos_input_feat = []
    pos_mask = []

    neg_rho_wrt_center = []
    neg_theta_wrt_center = []
    neg_input_feat = []
    neg_mask = []

    index = []
    pos_names = []
    neg_names = []

    idx_count = 0

    for ppi_pair_id in data_dirs:
        ppi_name = ppi_pair_id.rstrip('/').split('/')[-1]
        logger.debug('Loading %s', ppi_name)

        labels = np.load(os.path.join(ppi_pair_id, binder_name + '_sc_labels.npy'))

        # Take the median of the percentile 25 shape complementarity.
        mylabels = labels[0]
        labels = np.median(mylabels, axis=1)

        # Read the corresponding ply files.
        _, chains = ppi_name.split('_')
        # TODO make this work for different formats
        ply_fn1 = os.path.join(ply_dir, ppi_name + '_' + chains[:2] + '.ply')
        ply_fn2 = os.path.join(ply_dir, ppi_name + '_' + chains[2:] + '.ply')

        # pos_labels: points > max_sc_filt and >  min_sc_filt.
        # TODO add case for Nones
        pos_labels = np.where((labels < sc_max_filt) & (labels > sc_min_filt))[0]
        select = int(pos_surf_accept_probability * len(pos_labels))

        if select < 1:
            continue

        left = np.arange(len(pos_labels))
        np.random.shuffle(left)

        left = left[:select]
        left = pos_labels[left]

        v1 = pymesh.load_mesh(ply_fn1).vertices[left]
        v2 = pymesh.load_mesh(ply_fn2).vertices

        # For each point in v1, find the closest point in v2.
        kdt = cKDTree(v2)
        distances, right = kdt.query(v1)

        # Contact points: those within a cutoff distance.
        contact_points = np.where(distances < contact_distance)[0]
        k1 = left[contact_points]
        k2 = right[contact_points]

        # For negatives, get points in v2 far from p1.
        kdt = cKDTree(v1)
        dneg, _ = kdt.query(v2)
        k_neg2 = np.where(dneg > contact_distance)[0]

        assert len(k1) == len(k2)
        n_pos = len(k1)

        # Binder
        pid = binder_name
        for ii in k1:
            pos_names.append(f'{ppi_pair_id}_{pid}_{ii}')

        rho_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_rho_wrt_center.npy'))
        theta_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_theta_wrt_center.npy'))
        input_feat = np.load(os.path.join(ppi_pair_id, pid + '_input_feat.npy'))
        mask = np.load(os.path.join(ppi_pair_id, pid + '_mask.npy'))

        binder_rho_wrt_center.append(rho_wrt_center[k1])
        binder_theta_wrt_center.append(theta_wrt_center[k1])
        binder_input_feat.append(input_feat[k1])
        binder_mask.append(mask[k1])

        # Read pos
        pid = positive_name

        # Read as positives those points.
        rho_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_rho_wrt_center.npy'))
        theta_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_theta_wrt_center.npy'))
        input_feat = np.load(os.path.join(ppi_pair_id, pid+'_input_feat.npy'))
        mask = np.load(os.path.join(ppi_pair_id, pid+'_mask.npy'))

        pos_rho_wrt_center.append(rho_wrt_center[k2])
        pos_theta_wrt_center.append(theta_wrt_center[k2])
        pos_input_feat.append(input_feat[k2])
        pos_mask.append(mask[k2])

        # Get a set of negatives from  p2.
        np.random.shuffle(k_neg2)
        k_neg2 = k_neg2[:(len(k2))]

        assert len(k_neg2) == n_pos

        neg_rho_wrt_center.append(rho_wrt_center[k_neg2])
        neg_theta_wrt_center.append(theta_wrt_center[k_neg2])
        neg_input_feat.append(input_feat[k_neg2])
        neg_mask.append(mask[k_neg2])

        for ii in k_neg2:
            neg_names.append(f'{ppi_pair_id}_{pid}_{ii}')

        # Training, validation or test?
        index = np.append(index, np.arange(idx_count, idx_count + n_pos)).astype(int)
        idx_count += n_pos

    binder_rho_wrt_center = np.concatenate(binder_rho_wrt_center, axis=0)
    binder_theta_wrt_center = np.concatenate(binder_theta_wrt_center, axis=0)
    binder_input_feat = np.concatenate(binder_input_feat, axis=0)
    binder_mask = np.concatenate(binder_mask, axis=0)

    pos_rho_wrt_center = np.concatenate(pos_rho_wrt_center, axis=0)
    pos_theta_wrt_center = np.concatenate(pos_theta_wrt_center, axis=0)
    pos_input_feat = np.concatenate(pos_input_feat, axis=0)
    pos_mask = np.concatenate(pos_mask, axis=0)

    neg_rho_wrt_center = np.concatenate(neg_rho_wrt_center, axis=0)
    neg_theta_wrt_center = np.concatenate(neg_theta_wrt_center, axis=0)
    neg_input_feat = np.concatenate(neg_input_feat, axis=0)
    neg_mask = np.concatenate(neg_mask, axis=0)

    return (index,
            binder_rho_wrt_center, binder_theta_wrt_center, binder_input_feat, binder_mask,
            pos_rho_wrt_center, pos_theta_wrt_center, pos_input_feat, pos_mask,
            neg_rho_wrt_center, neg_theta_wrt_center, neg_input_feat, neg_mask)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    if args.seed is not None:
        logger.info('Seeding random state with seed %d', args.seed)
        np.random.seed(args.seed)

    logger.info('Aggregating training data')

    (
        training_idx,
        binder_rho_wrt_center, binder_theta_wrt_center, binder_input_feat, binder_mask,
        pos_rho_wrt_center, pos_theta_wrt_center, pos_input_feat, pos_mask,
        neg_rho_wrt_center, neg_theta_wrt_center, neg_input_feat, neg_mask
    ) = aggregate_data(args.training,
                       args.binder_name,
                       args.positive_name,
                       args.ply_dir,
                       args.contact_distance,
                       args.pos_surf_accept_probability,
                       args.sc_min_cutoff,
                       args.sc_max_cutoff)

    binder_input_feat = mask_input_feat(binder_input_feat, args.feat_mask)
    pos_input_feat = mask_input_feat(pos_input_feat, args.feat_mask)
    neg_input_feat = mask_input_feat(neg_input_feat, args.feat_mask)

    logger.debug('Number of training shapes: %d', len(training_idx))
    logger.debug('Read %d positive shapes', len(pos_rho_wrt_center))
    logger.debug('Read %d negative shapes', len(neg_rho_wrt_center))

    if args.validation_data:
        (
            val_idx,
            val_binder_rho_wrt_center, val_binder_theta_wrt_center, val_binder_input_feat, val_binder_mask,
            val_pos_rho_wrt_center, val_pos_theta_wrt_center, val_pos_input_feat, val_pos_mask,
            val_neg_rho_wrt_center, val_neg_theta_wrt_center, val_neg_input_feat, val_neg_mask,
        ) = aggregate_data(args.validation_data,
                           args.binder_name,
                           args.positive_name,
                           args.ply_dir,
                           args.contact_distance,
                           args.pos_surf_accept_probability,
                           args.sc_min_cutoff,
                           args.sc_max_cutoff)

        logger.debug('Number of validation shapes: %d', len(val_idx))

    learning_obj = MasifPPISearch(args.max_distance,
                                  n_thetas=16,
                                  n_rhos=5,
                                  n_rotations=16,
                                  idx_gpu='/gpu:0',
                                  feat_mask=args.feat_mask)

    if args.model:
        learning_obj.saver.restore(learning_obj.session, args.model)

    list_training_loss = []
    best_val_auc = 0

    pos_training_idx_copy = np.copy(training_idx)
    neg_training_idx_copy = np.copy(training_idx)

    logger.info('Number of iterations: %d', args.num_iterations)

    iter_pos_score = []
    iter_neg_score = []

    for num_iter in range(1, args.num_iterations + 1):
        logger.debug('Iterations number %d', num_iter)
        # Read dataset for training.
        np.random.shuffle(pos_training_idx_copy)
        np.random.shuffle(neg_training_idx_copy)

        c_pos_training_idx = pos_training_idx_copy[: args.batch_size // 4]
        c_neg_training_idx = neg_training_idx_copy[: args.batch_size // 4]

        # Features and theta are flipped for the binder in construct_batch (except for hydrophobicity).
        batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask = construct_batch(
            binder_rho_wrt_center,
            binder_theta_wrt_center,
            binder_input_feat,
            binder_mask,
            c_pos_training_idx,
            pos_rho_wrt_center,
            pos_theta_wrt_center,
            pos_input_feat,
            pos_mask,
            c_neg_training_idx,
            neg_rho_wrt_center,
            neg_theta_wrt_center,
            neg_input_feat,
            neg_mask,
        )

        assert len(batch_rho_coords) == args.batch_size
        assert len(batch_theta_coords) == args.batch_size
        assert len(batch_input_feat) == args.batch_size
        assert len(batch_mask) == args.batch_size

        feed_dict = {
            learning_obj.rho_coords: batch_rho_coords,
            learning_obj.theta_coords: batch_theta_coords,
            learning_obj.input_feat: batch_input_feat,
            learning_obj.mask: batch_mask,
            learning_obj.keep_prob: 0.5,
        }

        # Do not train during the first iteration
        if num_iter == 0:
            [score] = learning_obj.session.run(
                [learning_obj.score], feed_dict=feed_dict
            )
            training_loss = 0

        else:
            _, training_loss, _, score = learning_obj.session.run(
                [
                    learning_obj.optimizer,
                    learning_obj.data_loss,
                    learning_obj.norm_grad,
                    learning_obj.score,
                ],
                feed_dict=feed_dict,
            )

        n = len(score) // 2

        pos_score = score[:n]
        neg_score = score[n:]

        iter_pos_score = np.concatenate([pos_score, iter_pos_score], axis=0)
        iter_neg_score = np.concatenate([neg_score, iter_neg_score], axis=0)
        list_training_loss.append(training_loss)

        if num_iter % args.num_iter_eval == 0:
            logger.info('Evaluating at iteration %d', num_iter)

            logger.info('averaged training loss for the last %d iterations: %f',
                        args.num_iter_eval,
                        np.mean(list_training_loss))

            roc_auc = 1 - compute_roc_auc(iter_pos_score, iter_neg_score)
            logger.info('Training ROC-AUC: %f', roc_auc)

            logger.info('Mean training positive score: %f', np.mean(1.0 / iter_pos_score))
            logger.info('Mean training negative score: %f', np.mean(1.0 / iter_neg_score))

            iter_pos_score = []
            iter_neg_score = []
            list_training_loss = []

            if args.validation_data:
                pos_desc = compute_val_test_desc(
                    learning_obj,
                    val_idx,
                    val_pos_rho_wrt_center,
                    val_pos_theta_wrt_center,
                    val_pos_input_feat,
                    val_pos_mask,
                    batch_size=args.validation_batch_size,
                )

                binder_desc = compute_val_test_desc(
                    learning_obj,
                    val_idx,
                    val_binder_rho_wrt_center,
                    val_binder_theta_wrt_center,
                    val_binder_input_feat,
                    val_binder_mask,
                    batch_size=args.validation_batch_size,
                    flip=True,
                )

                neg_desc = compute_val_test_desc(
                    learning_obj,
                    val_idx,
                    val_neg_rho_wrt_center,
                    val_neg_theta_wrt_center,
                    val_neg_input_feat,
                    val_neg_mask,
                    batch_size=args.validation_batch_size,
                )

                neg_desc_2 = binder_desc.copy()

                # Simply shuffle negative descriptors.
                np.random.shuffle(neg_desc)

                # Compute val ROC AUC.
                pos_dists = compute_dists(pos_desc, binder_desc)
                neg_dists = compute_dists(neg_desc, neg_desc_2)

                val_auc = 1 - compute_roc_auc(pos_dists, neg_dists)

                logger.info('Validation ROC-AUC: %f', val_auc)

                logger.info('Mean validation positive score: %f', np.mean(pos_dists))
                logger.info('Mean validation negative score: %f', np.mean(neg_dists))

                if val_auc > best_val_auc:
                    logger.info('Lower validation ROC-AUC achieved, saving model...')

                    best_val_auc = val_auc
                    output_model = os.path.join(args.output, 'model')
                    learning_obj.saver.save(learning_obj.session, output_model)

    if args.validation_data is None:
        logger.info('Saving model...')
        output_model = os.path.join(args.output, 'model')
        learning_obj.saver.save(learning_obj.session, output_model)


if __name__ == '__main__':
    main()
