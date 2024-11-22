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
'''Functions for creating surface representations and adding features to the surfaces.'''
import numpy as np
import trimesh

from scipy.spatial import cKDTree  # pylint: disable = no-name-in-module


def extract_patch_and_coord(vix, shape, coord, max_distance, max_vertices, patch_indices=False):
    '''From a full shape in a full protein, extract a patch around a vertex.

    If patch_indices = True, then store the indices of all neighbors.

    '''
    # Member vertices are nonzero elements
    _, j = coord[np.int(vix), : coord.shape[1] // 2].nonzero()

    d = np.squeeze(np.asarray(coord[np.int(vix), : coord.shape[1] // 2].todense()))
    j = np.where((d < max_distance) & (d > 0))[0]

    max_dist_tmp = max_distance

    while len(j) > max_vertices:
        max_dist_tmp = max_dist_tmp * 0.95
        j = np.where((d < max_dist_tmp) & (d > 0))[0]

    d = d[j]

    patch = {}
    patch['X'] = shape['X'][0][j]
    patch['Y'] = shape['Y'][0][j]
    patch['Z'] = shape['Z'][0][j]
    patch['charge'] = shape['charge'][0][j]
    patch['hbond'] = shape['hbond'][0][j]
    patch['normal'] = shape['normal'][:, j]
    patch['shape_index'] = shape['shape_index'][0][j]

    if 'hphob' in shape:
        patch['hphob'] = shape['hphob'][0][j]

    patch['center'] = np.argmin(d)

    j_theta = j + coord.shape[1] // 2
    theta = np.squeeze(np.asarray(coord[np.int(vix), j_theta].todense()))
    coord = np.concatenate([d, theta], axis=0)

    if patch_indices:
        return patch, coord, j

    else:
        return patch, coord


def compute_shape_complementarity(ply_fn1, ply_fn2,
                                  neigh1, neigh2,
                                  rho1, rho2,
                                  mask1, mask2,
                                  sc_w,
                                  sc_interaction_cutoff,
                                  sc_radius):
    '''Compute the shape complementarity between all pairs of patches.

    neigh1 and neigh2 are the precomputed indices; rho1 and rho2 their distances.

    Args:
        ply_fnX: path to the ply file of the surface of protein X=1 and X=2
        neighX, rhoX, maskX: (N,max_vertices_per_patch) matrices with the indices of the neighbors, the distances to the
                             center and the mask

        Returns:
            vX_sc (2,N,10) matrix with the shape complementarity (shape complementarity 25 and 50) of each vertex to its
            nearest neighbor in the other protein, in 10 rings.

    '''
    mesh1 = trimesh.load(ply_fn1)
    n1 = mesh1.vertex_normals

    mesh2 = trimesh.load(ply_fn2)
    n2 = mesh2.vertex_normals

    w = sc_w
    int_cutoff = sc_interaction_cutoff
    radius = sc_radius
    num_rings = 10

    scales = np.arange(0, radius, radius/10)
    scales = np.append(scales, radius)

    v1 = mesh1.vertices
    v2 = mesh2.vertices

    v1_sc = np.zeros((2, len(v1), 10))
    v2_sc = np.zeros((2, len(v2), 10))

    # Find all interface vertices
    kdt = cKDTree(v2)
    d, nearest_neighbors_v1_to_v2 = kdt.query(v1)

    # Interface vertices in v1
    interface_vertices_v1 = np.where(d < int_cutoff)[0]

    # Go through every interface vertex.
    for cv1_iiix in range(len(interface_vertices_v1)):
        cv1_ix = interface_vertices_v1[cv1_iiix]

        assert d[cv1_ix] < int_cutoff

        # First shape complementarity s1->s2 for the entire patch
        patch_idxs1 = np.where(mask1[cv1_ix] == 1)[0]
        neigh_cv1 = np.array(neigh1[cv1_ix])[patch_idxs1]

        # Find the point cv2_ix in s2 that is closest to cv1_ix
        cv2_ix = nearest_neighbors_v1_to_v2[cv1_ix]

        patch_idxs2 = np.where(mask2[cv2_ix] == 1)[0]
        neigh_cv2 = np.array(neigh2[cv2_ix])[patch_idxs2]

        patch_v1 = v1[neigh_cv1]
        patch_v2 = v2[neigh_cv2]

        patch_n1 = n1[neigh_cv1]
        patch_n2 = n2[neigh_cv2]

        patch_kdt = cKDTree(patch_v1)
        p_dists_v2_to_v1, p_nearest_neighbor_v2_to_v1 = patch_kdt.query(patch_v2)

        patch_kdt = cKDTree(patch_v2)
        p_dists_v1_to_v2, p_nearest_neighbor_v1_to_v2 = patch_kdt.query(patch_v1)

        # First v1->v2
        neigh_cv1_p = p_nearest_neighbor_v1_to_v2

        comp1 = [np.dot(patch_n1[x], -patch_n2[neigh_cv1_p][x]) for x in range(len(patch_n1))]
        comp1 = np.multiply(comp1, np.exp(-w * np.square(p_dists_v1_to_v2)))

        # Use 10 rings such that each ring has equal weight in shape complementarity
        comp_rings1_25 = np.zeros(num_rings)
        comp_rings1_50 = np.zeros(num_rings)

        patch_rho1 = np.array(rho1[cv1_ix])[patch_idxs1]
        for ring in range(num_rings):
            # scale = scales[ring]
            members = np.where((patch_rho1 >= scales[ring]) & (patch_rho1 < scales[ring + 1]))

            if len(members[0]) == 0:
                comp_rings1_25[ring] = 0.0
                comp_rings1_50[ring] = 0.0

            else:
                comp_rings1_25[ring] = np.percentile(comp1[members], 25)
                comp_rings1_50[ring] = np.percentile(comp1[members], 50)

        # Now v2->v1
        neigh_cv2_p = p_nearest_neighbor_v2_to_v1

        comp2 = [np.dot(patch_n2[x], -patch_n1[neigh_cv2_p][x]) for x in range(len(patch_n2))]
        comp2 = np.multiply(comp2, np.exp(-w * np.square(p_dists_v2_to_v1)))

        # Use 10 rings such that each ring has equal weight in shape complementarity
        comp_rings2_25 = np.zeros(num_rings)
        comp_rings2_50 = np.zeros(num_rings)

        # Apply mask to patch rho coordinates.
        patch_rho2 = np.array(rho2[cv2_ix])[patch_idxs2]

        for ring in range(num_rings):
            # scale = scales[ring]
            members = np.where((patch_rho2 >= scales[ring]) & (patch_rho2 < scales[ring + 1]))

            if len(members[0]) == 0:
                comp_rings2_25[ring] = 0.0
                comp_rings2_50[ring] = 0.0

            else:
                comp_rings2_25[ring] = np.percentile(comp2[members], 25)
                comp_rings2_50[ring] = np.percentile(comp2[members], 50)

        v1_sc[0, cv1_ix, :] = comp_rings1_25
        v2_sc[0, cv2_ix, :] = comp_rings2_25
        v1_sc[1, cv1_ix, :] = comp_rings1_50
        v2_sc[1, cv2_ix, :] = comp_rings2_50

    return v1_sc, v2_sc


def normalize_electrostatics(in_elec):
    '''Normalize electrostatics to a value between -1 and 1.'''
    elec = np.copy(in_elec)

    upper_threshold = 3
    lower_threshold = -3

    elec[elec > upper_threshold] = upper_threshold
    elec[elec < lower_threshold] = lower_threshold

    elec = elec - lower_threshold
    elec = elec / (upper_threshold - lower_threshold)

    elec = 2 * elec - 1

    return elec


def mean_normal_center_patch(d, n, r):
    '''Function to compute the mean normal of vertices within r radius of the center of the patch.'''
    c_normal = [n[i] for i in range(len(d)) if d[i] <= r]

    mean_normal = np.mean(c_normal, axis=0, keepdims=True).T
    mean_normal = mean_normal / np.linalg.norm(mean_normal)

    return np.squeeze(mean_normal)


def compute_ddc(patch_v, patch_n, patch_cp, patch_rho):
    '''Compute the distance dependent curvature, Yin et al PNAS 2009

    Args:
        patch_v: the patch vertices
        patch_n: the patch normals
        patch_cp: the index of the central point of the patch
        patch_rho: the geodesic distance to all members.

    Returns:
        a vector with the ddc for each point in the patch

    '''
    n = patch_n
    r = patch_v
    i = patch_cp

    # Compute the mean normal 2.5A around the center point
    ni = mean_normal_center_patch(patch_rho, n, 2.5)
    dij = np.linalg.norm(r - r[i], axis=1)

    # Compute the step function sf:
    sf = r + n
    sf = sf - (ni + r[i])
    sf = np.linalg.norm(sf, axis=1)
    sf = sf - dij

    sf[sf > 0] = 1
    sf[sf < 0] = -1
    sf[sf == 0] = 0

    # Compute the curvature between i and j
    dij[dij == 0] = 1e-8

    kij = np.divide(np.linalg.norm(n - ni, axis=1), dij)
    kij = np.multiply(sf, kij)

    # Ignore any values greater than 0.7 and any values smaller than 0.7
    kij[kij > 0.7] = 0
    kij[kij < -0.7] = 0

    return kij
