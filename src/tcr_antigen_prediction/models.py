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
'''Models used in package.'''
import logging

import torch
import torch.nn as nn
import numpy as np

logger = logging.getLogger(__name__)


class MasifPPISearch(nn.Module):
    '''The neural network model to classify two patches into binders or not binders.

    Note: Features for binder should be flipped before feeding to the NN.

    '''

    def __init__(self,
                 max_rho: float,
                 n_thetas: int = 16,
                 n_rhos: int = 5,
                 n_rotations: int = 16,
                 feat_mask: list[int] | None = None):
        super().__init__()

        if feat_mask is None:
            feat_mask = [1.0, 1.0, 1.0, 1.0, 1.0]

        self.feat_mask = feat_mask

        self.n_rhos = n_rhos
        self.n_thetas = n_thetas
        self.n_feats = int(sum(feat_mask))

        initial_coords = compute_initial_coordinates(max_rho, n_rhos, n_thetas)

        mu_rho_initial = torch.from_numpy(np.expand_dims(initial_coords[:, 0], 0).astype('float32'))
        mu_theta_initial = torch.from_numpy(np.expand_dims(initial_coords[:, 1], 0).astype('float32'))

        self.convs = nn.ParameterList([PolarConvolutionalLayer(mu_rho_initial,
                                                               mu_theta_initial,
                                                               max_rho,
                                                               n_rhos,
                                                               n_thetas,
                                                               n_rotations) for _ in range(self.n_feats)])
        self.output_fc = nn.Linear(n_thetas * n_rhos * self.n_feats, n_thetas * n_rhos)

        # # compute data loss
        # self.n_patches = tf.shape(self.global_desc)[0] // 4
        # self.data_loss = self.compute_data_loss()

        # # definition of the solver
        # self.optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(self.data_loss)
        # self.var_grad = tf.gradients(self.data_loss, tf.trainable_variables())

        # for k in range(len(self.var_grad)):
        #     if self.var_grad[k] is None:
        #         logger.debug(tf.trainable_variables()[k])

        # self.norm_grad = self.frobenius_norm(tf.concat([tf.reshape(g, [-1]) for g in self.var_grad], 0))

        # # Create a session for running Ops on the Graph.
        # config = tf.ConfigProto(allow_soft_placement=True)
        # config.gpu_options.allow_growth = True  # pylint: disable = no-member

        # self.session = tf.Session(config=config)
        # self.saver = tf.train.Saver()

        # # Run the Op to initialize the variables.
        # init = tf.global_variables_initializer()

        # self.session.run(init)
        # self.count_number_parameters()

    def forward(self,
                input_feats: torch.Tensor,
                rho_coords: torch.Tensor,
                theta_coords: torch.Tensor,
                mask: torch.Tensor):
        feat_descs = []
        for i in range(self.n_feats):
            logger.debug('Feature %d', i)

            input_i = torch.unsqueeze(input_feats[:, :, i], 2)

            desc = self.convs[i](input_i,
                                 rho_coords,
                                 theta_coords,
                                 mask)

            feat_descs.append(desc)

        descs = torch.stack(feat_descs, dim=1)
        descs = torch.reshape(descs, [-1, self.n_rhos * self.n_thetas * self.n_feats])
        descs = self.output_fc(descs)

        return descs

    # def frobenius_norm(self, tensor):
    #     square_tensor = tf.square(tensor)
    #     tensor_sum = tf.reduce_sum(square_tensor)
    #     frobenius_norm = tf.sqrt(tensor_sum)

    #     return frobenius_norm

    # def build_sparse_matrix_softmax(self, idx_non_zero_values, mat, dense_shape_out):
    #     out_mat = tf.SparseTensorValue(idx_non_zero_values, tf.squeeze(mat), dense_shape_out)
    #     out_mat = tf.sparse_reorder(out_mat)  # n_edges x n_edges
    #     out_mat = tf.sparse_softmax(out_mat)

    #     return out_mat

    # def compute_data_loss_cross_entropy(self, pos, neg):
    #     '''Softmax cross entropy'''
    #     epsilon = tf.constant(value=0.00001)
    #     logit = tf.nn.softmax([pos, neg])

    #     self.softmax_debug = logit

    #     cross_entropy = -(tf.log(logit[1] + epsilon) - tf.log(logit[0] + epsilon))
    #     return cross_entropy

    # def compute_data_loss(self, pos_thresh=0.0, neg_thresh=10):
    #     '''Data loss. Values above 10 are ignored.'''
    #     self.global_desc_pos = tf.gather(self.global_desc, tf.range(0, self.n_patches))
    #     self.global_desc_binder = tf.gather(self.global_desc, tf.range(self.n_patches, 2 * self.n_patches))
    #     self.global_desc_neg = tf.gather(self.global_desc, tf.range(2 * self.n_patches, 3 * self.n_patches))
    #     self.global_desc_neg_2 = tf.gather(self.global_desc, tf.range(3 * self.n_patches, 4 * self.n_patches))

    #     pos_distances = tf.reduce_sum(tf.square(self.global_desc_binder - self.global_desc_pos), 1)
    #     neg_distances = tf.reduce_sum(tf.square(self.global_desc_neg - self.global_desc_neg_2), 1)

    #     self.score = tf.concat([pos_distances, neg_distances], 0)

    #     pos_distances = tf.nn.relu(tf.reduce_sum(tf.square(self.global_desc_binder - self.global_desc_pos), 1)
    #                                - pos_thresh)
    #     neg_distances = tf.nn.relu(-tf.reduce_sum(tf.square(self.global_desc_neg - self.global_desc_neg_2), 1)
    #                                + neg_thresh)

    #     pos_mean, pos_std = tf.nn.moments(pos_distances, [0])
    #     neg_mean, neg_std = tf.nn.moments(neg_distances, [0])

    #     data_loss = pos_std + neg_std + pos_mean + neg_mean

    #     return data_loss


def compute_initial_coordinates(max_rho, n_rhos, n_thetas):
    range_rho = [0.0, max_rho]
    range_theta = [0, 2 * np.pi]

    grid_rho = np.linspace(range_rho[0], range_rho[1], num=n_rhos + 1)
    grid_rho = grid_rho[1:]

    grid_theta = np.linspace(range_theta[0], range_theta[1], num=n_thetas + 1)
    grid_theta = grid_theta[:-1]

    grid_rho_, grid_theta_ = np.meshgrid(grid_rho, grid_theta, sparse=False)

    # the traspose here is needed to have the same behaviour as Matlab code
    grid_rho_ = grid_rho_.T
    grid_theta_ = grid_theta_.T

    grid_rho_ = grid_rho_.flatten()
    grid_theta_ = grid_theta_.flatten()

    coords = np.concatenate((grid_rho_[None, :], grid_theta_[None, :]), axis=0)
    coords = coords.T  # every row contains the coordinates of a grid intersection

    logger.debug(coords.shape)

    return coords


class PolarConvolutionalLayer(nn.Module):
    def __init__(self,
                 mu_rho_initial: np.ndarray,
                 mu_theta_initial: np.ndarray,
                 max_rho: float,
                 n_rhos: int,
                 n_thetas: int,
                 n_rotations: int,
                 mean_gauss_activation: bool = True,
                 eps: float = 1e-5):
        super().__init__()

        self.n_rotations = n_rotations

        self.gauss_activation = GaussianActivation(mu_rho_initial,
                                                   mu_theta_initial,
                                                   max_rho,
                                                   n_rhos,
                                                   n_thetas,
                                                   mean_gauss_activation,
                                                   eps)

        self.weight = nn.Parameter(torch.empty(n_thetas * n_rhos, n_thetas * n_rhos))
        nn.init.xavier_uniform_(self.weight)

        self.bias = nn.Parameter(torch.zeros([n_thetas * n_rhos]))

    def forward(self,
                input_feat: torch.Tensor,
                rho_coords: torch.Tensor,
                theta_coords: torch.Tensor,
                mask: torch.Tensor):
        all_conv_feat = []
        for k in range(self.n_rotations):
            logger.debug('Rotation %d', k)

            rho_coords = torch.reshape(rho_coords, [-1, 1])       # batch_size * n_vertices
            theta_coords = torch.reshape(theta_coords, [-1, 1])  # batch_size * n_vertices

            theta_coords = torch.add(theta_coords, k * 2 * np.pi / self.n_rotations)
            theta_coords = torch.remainder(theta_coords, 2 * np.pi)

            gauss_desc = self.gauss_activation(input_feat, rho_coords, theta_coords, mask)

            conv_feat = torch.matmul(gauss_desc, self.weight) + self.bias    # batch_size, 80

            all_conv_feat.append(conv_feat)

        all_conv_feat = torch.stack(all_conv_feat)

        conv_feat = torch.max(all_conv_feat, dim=0).values
        conv_feat = torch.relu(conv_feat)

        return conv_feat


class GaussianActivation(nn.Module):
    def __init__(self,
                 mu_rho_initial: torch.Tensor,
                 mu_theta_initial: torch.Tensor,
                 max_rho: float,
                 n_rhos: int,
                 n_thetas: int,
                 mean_gauss_activation: bool = True,
                 eps: float = 1e-5):
        super().__init__()

        self.n_rhos = n_rhos
        self.n_thetas = n_thetas

        self.mean_gauss_activation = mean_gauss_activation
        self.eps = eps

        # order of the spectral filters
        sigma_rho_init = max_rho / 8  # in MoNet was 0.005 with max radius = 0.04 (i.e. 8 times smaller)
        sigma_theta_init = 1.0

        self.mu_rho = nn.Parameter(mu_rho_initial)      # 1, n_gauss
        self.mu_theta = nn.Parameter(mu_theta_initial)  # 1, n_gauss

        self.sigma_rho = nn.Parameter(torch.ones_like(mu_rho_initial) * sigma_rho_init)            # 1, n_gauss
        self.sigma_theta = nn.Parameter(torch.ones_like(mu_theta_initial) * sigma_theta_init)      # 1, n_gauss

    def forward(self,
                input_feat: torch.Tensor,
                rho_coords: torch.Tensor,
                theta_coords: torch.Tensor,
                mask: torch.Tensor):
        n_samples = mask.shape[0]
        n_vertices = mask.shape[1]

        rho_coords = torch.exp(-torch.square(rho_coords - self.mu_rho) / (torch.square(self.sigma_rho) + self.eps))
        theta_coords = torch.exp(-torch.square(theta_coords - self.mu_theta)
                                 / (torch.square(self.sigma_theta) + self.eps))

        gauss_activations = torch.mul(rho_coords, theta_coords)       # batch_size * n_vertices, n_gauss
        gauss_activations = torch.reshape(gauss_activations,          # batch_size, n_vertices, n_gauss
                                          [n_samples, n_vertices, -1])
        gauss_activations = torch.multiply(gauss_activations, mask)

        if self.mean_gauss_activation:
            # batch_size, n_vertices, n_gauss
            gauss_activations = torch.div(gauss_activations,
                                          torch.sum(gauss_activations, dim=1, keepdim=True) + self.eps)

        gauss_activations = torch.unsqueeze(gauss_activations, 2)    # batch_size, n_vertices, 1, n_gauss,

        input_feat_ = torch.unsqueeze(input_feat, 3)                 # batch_size, n_vertices, n_feats, 1

        gauss_desc = torch.multiply(gauss_activations, input_feat_)  # batch_size, n_vertices, n_feats, n_gauss,
        gauss_desc = torch.sum(gauss_desc, 1)                        # batch_size, n_feats, n_gauss,
        gauss_desc = torch.reshape(gauss_desc,                       # batch_size, 80
                                   [n_samples, self.n_thetas * self.n_rhos])

        return gauss_desc
