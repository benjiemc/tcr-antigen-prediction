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
'''Command line application to compute numerical representations of binding proteins.'''
import argparse
import logging
import os
import sys

import numpy as np
import torch

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.models import MasifPPISearch

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--model', required=True, help='path to trained model')
parser.add_argument('--output', '-o', required=True, help='path to output')
parser.add_argument('input', nargs='+', help='path to input directories to compute descriptors')

add_logging_arguments(parser)

np.random.seed(0)

params = {
    'max_shape_size': 200,
    'max_distance': 12.0,
    'feat_mask': [1.0] * 5,
    'max_sc_filt': 1.0,
    'min_sc_filt': 0.5,
    'pos_surf_accept_probability': 1.0,
    'pos_interface_cutoff': 1.0,
    'range_val_samples': 0.9,
    'sc_radius': 12.0,
    'sc_interaction_cutoff': 1.5,
    'sc_w': 0.25,
}


def mask_input_feat(input_feat: np.ndarray, mask: np.ndarray):
    '''Apply mask to input_feat'''
    mymask = np.where(mask == 0.0)[0]
    return np.delete(input_feat, mymask, axis=2)


def construct_batch(c_idx: np.ndarray,
                    input_feat: np.ndarray,
                    rho_wrt_center: np.ndarray,
                    theta_wrt_center: np.ndarray,
                    mask: np.ndarray,
                    flip=False):
    '''Construct a batch from the indices.'''
    batch_input_feat = input_feat[c_idx]
    batch_rho_coords = np.expand_dims(rho_wrt_center[c_idx], 2)
    batch_theta_coords = np.expand_dims(theta_wrt_center[c_idx], 2)
    batch_mask = np.expand_dims(mask[c_idx], 2)

    # Flip features and theta (except hydrophobicity)
    if flip:
        batch_input_feat = -batch_input_feat
        batch_theta_coords = 2 * np.pi - batch_theta_coords

        assert len(batch_input_feat.shape) == 3

        # Hydrophobicity is not flipped. -- Fix this.
        if batch_input_feat.shape[2] == 5 or batch_input_feat.shape[2] == 3:
            batch_input_feat[:, :, -1] = -batch_input_feat[:, :, -1]

    return batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask


def compute_descriptors(model: torch.nn.Module,
                        idx: np.ndarray,
                        input_feat: np.ndarray,
                        rho_wrt_center: np.ndarray,
                        theta_wrt_center: np.ndarray,
                        mask: np.ndarray,
                        batch_size: int = 100,
                        flip: bool = False,
                        device: str = 'cpu'):
    all_descs = []
    num_batches = int(np.ceil(float(len(idx)) / float(batch_size)))

    # Compute all desc for positive shapes.
    for kk in range(num_batches):
        c_idx = idx[np.arange(kk * batch_size, min((kk + 1) * batch_size, len(idx)))]

        batch = construct_batch(c_idx, input_feat, rho_wrt_center, theta_wrt_center, mask, flip=flip)
        batch = [torch.from_numpy(item).to(device) for item in batch]

        desc = model(*batch)
        desc = np.squeeze(desc.numpy())

        if len(desc.shape) == 1:
            desc = np.expand_dims(desc, 0)

        all_descs.append(desc)

    if len(all_descs) > 1:
        all_descs = np.concatenate(all_descs, axis=0)

    else:
        all_descs = all_descs[0]

    return all_descs


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    model = MasifPPISearch(params['max_distance'],
                           n_thetas=16,
                           n_rhos=5,
                           n_rotations=16,
                           feat_mask=params['feat_mask'])
    model.load_state_dict(torch.load(args.model, weights_only=True))

    if torch.cuda.is_available():
        device = torch.device('cuda')

    else:
        logger.warning('CUDA is not avaible, running on CPU which can be slow')
        device = torch.device('cpu')

    model.to(device)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    num_pairs = len(args.input)

    for num, ppi_pair_id in enumerate(args.input, 1):
        ppi_pair_id_name = os.path.basename(ppi_pair_id)

        logger.info('Working on: %s - %d of %d', ppi_pair_id_name, num, num_pairs)

        out_desc_dir = os.path.join(args.output, ppi_pair_id_name)

        if not os.path.exists(out_desc_dir):
            os.mkdir(out_desc_dir)

        for pid in 'p1', 'p2':
            input_feat = np.load(os.path.join(ppi_pair_id, pid + '_input_feat.npy'))
            input_feat = mask_input_feat(input_feat, np.array(params['feat_mask']))

            rho_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_rho_wrt_center.npy'))
            theta_wrt_center = np.load(os.path.join(ppi_pair_id, pid + '_theta_wrt_center.npy'))

            mask = np.load(os.path.join(ppi_pair_id, pid + '_mask.npy'))

            idx = np.array(range(len(rho_wrt_center)))

            desc_str = compute_descriptors(model,
                                           idx,
                                           input_feat,
                                           rho_wrt_center,
                                           theta_wrt_center,
                                           mask,
                                           batch_size=1000,
                                           flip=False,
                                           device=device)

            desc_flip = compute_descriptors(model,
                                            idx,
                                            input_feat,
                                            rho_wrt_center,
                                            theta_wrt_center,
                                            mask,
                                            batch_size=1000,
                                            flip=True,
                                            device=device)

            np.save(os.path.join(out_desc_dir, f'{pid}_desc_straight.npy'), desc_str)
            np.save(os.path.join(out_desc_dir, f'{pid}_desc_flipped.npy'), desc_flip)


if __name__ == '__main__':
    main()
