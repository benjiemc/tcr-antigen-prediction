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
'''Functions for creating meshes and adding features to those meshes from PDB structures.'''
import os
import logging
import tempfile
from subprocess import Popen, PIPE

import numpy as np
import pymesh
from Bio.PDB import PDBParser, Selection, NeighborSearch
from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral
from numpy.matlib import repmat
from sklearn.neighbors import KDTree

from tcr_antigen_prediction.io import output_pdb_as_xyzrn, read_msms
from tcr_antigen_prediction.chemistry import (DONOR_ATOM,
                                              POLAR_HYDROGENS,
                                              ACCEPTOR_ANGLE_ATOM,
                                              ACCEPTOR_PLANE_ATOM,
                                              HBOND_STD_DEV)

logger = logging.getLogger(__name__)

EPSILON = 1.0e-6


def compute_charges(pdb_filename, vertices, names):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_filename, pdb_filename)

    residues = {}
    for res in struct.get_residues():
        chain_id = res.get_parent().get_id()

        if chain_id == '':
            chain_id = ' '

        residues[(chain_id, res.get_id())] = res

    atoms = Selection.unfold_entities(struct, 'A')
    satisfied_co, satisfied_hn = compute_statisfied_co_hn(atoms)

    charge = np.array([0.0] * len(vertices))

    # Go over every vertex
    for ix, name in enumerate(names):
        fields = name.split('_')
        chain_id = fields[0]

        if chain_id == '':
            chain_id = ' '

        if fields[2] == 'x':
            fields[2] = ' '

        res_id = (' ', int(fields[1]), fields[2])
        atom_name = fields[4]

        # Ignore atom if it is BB and it is already satisfied.
        if atom_name == 'H' and res_id in satisfied_hn:
            continue

        if atom_name == 'O' and res_id in satisfied_co:
            continue

        # Compute the charge of the vertex
        charge[ix] = _compute_charge_helper(atom_name, residues[(chain_id, res_id)], vertices[ix])

    return charge


def _compute_charge_helper(atom_name, res, v):
    '''Compute the charge of a vertex in a residue.'''
    # Check if it is a polar hydrogen.
    if is_polar_hydrogen(atom_name, res):
        donor_atom_name = DONOR_ATOM[atom_name]

        a = res[donor_atom_name].get_coord()  # N/O
        b = res[atom_name].get_coord()  # H

        # Donor-H is always 180.0 degrees, = pi
        angle_deviation = compute_angle_deviation(a, b, v, np.pi)
        angle_penalty = compute_angle_penalty(angle_deviation)

        return 1.0 * angle_penalty

    # Check if it is an acceptor oxygen or nitrogen
    elif is_acceptor_atom(atom_name, res):
        acceptor_atom = res[atom_name]

        a = res[ACCEPTOR_ANGLE_ATOM[atom_name]].get_coord()
        b = acceptor_atom.get_coord()

        # TODO: This should not be 120 for all atoms, i.e. for HIS it should be ~125.0
        angle_deviation = compute_angle_deviation(a, b, v, 2 * np.pi / 3)  # 120 degress for acceptor
        angle_penalty = compute_angle_penalty(angle_deviation)

        plane_penalty = 1.0

        if atom_name in ACCEPTOR_PLANE_ATOM:
            d = res[ACCEPTOR_PLANE_ATOM[atom_name]].get_coord()

            plane_deviation = compute_plane_deviation(d, a, b, v)
            plane_penalty = compute_angle_penalty(plane_deviation)

        return -1.0 * angle_penalty * plane_penalty

    return 0.0


def compute_angle_deviation(a, b, c, theta):
    '''Compute the absolute value of the deviation from theta'''
    return abs(calc_angle(Vector(a), Vector(b), Vector(c)) - theta)


def compute_plane_deviation(a, b, c, d):
    '''Compute the angle deviation from a plane'''
    dih = calc_dihedral(Vector(a), Vector(b), Vector(c), Vector(d))

    dev1 = abs(dih)
    dev2 = np.pi - abs(dih)

    return min(dev1, dev2)


def compute_angle_penalty(angle_deviation):
    '''Angle_deviation from ideal value.

    TODO: do a more data-based solution

    '''
    return max(0.0, 1.0 - (angle_deviation / (HBOND_STD_DEV)) ** 2)  # Standard deviation: HBOND_STD_DEV


def is_polar_hydrogen(atom_name, res):
    return atom_name in POLAR_HYDROGENS[res.get_resname()]


def is_acceptor_atom(atom_name, res):
    if atom_name.startswith('O'):
        return True

    else:
        if res.get_resname() == 'HIS':
            if atom_name == 'ND1' and 'HD1' not in res:
                return True

            if atom_name == 'NE2' and 'HE2' not in res:
                return True

    return False


def compute_statisfied_co_hn(atoms):
    '''Compute the list of backbone C=O:H-N that are satisfied. These will be ignored.'''
    ns = NeighborSearch(atoms)

    satisfied_co = set()
    satisfied_hn = set()

    for atom1 in atoms:
        res1 = atom1.get_parent()

        if atom1.get_id() == 'O':
            neigh_atoms = ns.search(atom1.get_coord(), 2.5, level='A')

            for atom2 in neigh_atoms:
                if atom2.get_id() == 'H':
                    res2 = atom2.get_parent()

                    # Ensure they belong to different residues.
                    if res2.get_id() != res1.get_id():
                        # Compute the angle N-H:O, ideal value is 180 (but in
                        # helices it is typically 160) 180 +-30 = pi
                        angle_n_h_o_dev = compute_angle_deviation(res2['N'].get_coord(),
                                                                  atom2.get_coord(),
                                                                  atom1.get_coord(),
                                                                  np.pi)

                        # Compute angle H:O=C, ideal value is ~160 +- 20 = 8*pi/9
                        angle_h_o_c_dev = compute_angle_deviation(atom2.get_coord(),
                                                                  atom1.get_coord(),
                                                                  res1['C'].get_coord(),
                                                                  8 * np.pi / 9)

                        # Allowed deviations: 30 degrees (pi/6) and 20 degrees (pi/9)
                        if (angle_n_h_o_dev - np.pi / 6 < 0) and (angle_h_o_c_dev - np.pi / 9 < 0.0):
                            satisfied_co.add(res1.get_id())
                            satisfied_hn.add(res2.get_id())

    return satisfied_co, satisfied_hn


def assign_charges_to_new_mesh(new_vertices, old_vertices, old_charges):
    '''Compute the charge of a new mesh, based on the charge of an old mesh.

    Use the top vertex in distance, for now (later this should be smoothed over 3 or 4 vertices)

    '''
    dataset = old_vertices
    testset = new_vertices

    new_charges = np.zeros(len(new_vertices))

    num_inter = 4  # Number of interpolation features

    # Assign k old vertices to each new vertex.
    kdt = KDTree(dataset)
    dists, result = kdt.query(testset, k=num_inter)

    # Square the distances (as in the original pyflann)
    dists = np.square(dists)

    # The size of result is the same as new_vertices
    for vi_new in range(len(result)):
        vi_old = result[vi_new]
        dist_old = dists[vi_new]

        # If one vertex is right on top, ignore the rest.
        if dist_old[0] == 0.0:
            new_charges[vi_new] = old_charges[vi_old[0]]
            continue

        total_dist = np.sum(1 / dist_old)
        for i in range(num_inter):
            new_charges[vi_new] += old_charges[vi_old[i]] * (1 / dist_old[i]) / total_dist

    return new_charges


def compute_msms(pdb_file,  protonate=True):
    with tempfile.TemporaryDirectory() as temp_dir:
        base_name = os.path.join(temp_dir, os.path.basename(pdb_file).rsplit('.', 1)[0])
        out_xyzrn = base_name + '.xyzrn'

        if protonate:
            output_pdb_as_xyzrn(pdb_file, out_xyzrn)

        # Now run MSMS on xyzrn file
        args = ['msms',
                '-density', '3.0',
                '-hdensity', '3.0',
                '-probe', '1.5',
                '-if', out_xyzrn,
                '-of', base_name,
                '-af', base_name]

        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        _, _ = p2.communicate()

        vertices, faces, normals, names = read_msms(base_name)
        areas = {}
        with open(base_name + '.area') as ses_file:
            next(ses_file)  # ignore header line

            for line in ses_file:
                fields = line.split()
                areas[fields[3]] = fields[1]

    return vertices, faces, normals, names, areas


kd_scale = {
    'ILE': 4.5,
    'VAL': 4.2,
    'LEU': 3.8,
    'PHE': 2.8,
    'CYS': 2.5,
    'MET': 1.9,
    'ALA': 1.8,
    'GLY': -0.4,
    'THR': -0.7,
    'SER': -0.8,
    'TRP': -0.9,
    'TYR': -1.3,
    'PRO': -1.6,
    'HIS': -3.2,
    'GLU': -3.5,
    'GLN': -3.5,
    'ASP': -3.5,
    'ASN': -3.5,
    'LYS': -3.9,
    'ARG': -4.5,
}
'''Kyte Doolittle scale'''


def compute_hydrophobicity(names):
    '''For each vertex in names, compute'''
    hp = np.zeros(len(names))
    for ix, name in enumerate(names):
        aa = name.split('_')[3]
        hp[ix] = kd_scale[aa]

    return hp


def fix_mesh(mesh, resolution=1.0, detail='normal'):
    bbox_min, bbox_max = mesh.bbox
    diag_len = np.linalg.norm(bbox_max - bbox_min)

    if detail == 'normal':
        target_len = diag_len * 5e-3

    elif detail == 'high':
        target_len = diag_len * 2.5e-3

    elif detail == 'low':
        target_len = diag_len * 1e-2

    target_len = resolution

    mesh, _ = pymesh.remove_duplicated_vertices(mesh, 0.001)

    count = 0

    logger.info('Removing degenerated triangles')
    mesh, _ = pymesh.remove_degenerated_triangles(mesh, 100)
    mesh, _ = pymesh.split_long_edges(mesh, target_len)
    num_vertices = mesh.num_vertices

    while True:
        mesh, _ = pymesh.collapse_short_edges(mesh, 1e-6)
        mesh, _ = pymesh.collapse_short_edges(mesh, target_len, preserve_feature=True)
        mesh, _ = pymesh.remove_obtuse_triangles(mesh, 150.0, 100)

        if mesh.num_vertices == num_vertices:
            break

        num_vertices = mesh.num_vertices
        count += 1

        if count > 10:
            break

    mesh = pymesh.resolve_self_intersection(mesh)
    mesh, _ = pymesh.remove_duplicated_faces(mesh)
    mesh = pymesh.compute_outer_hull(mesh)
    mesh, _ = pymesh.remove_duplicated_faces(mesh)
    mesh, _ = pymesh.remove_obtuse_triangles(mesh, 179.0, 5)
    mesh, _ = pymesh.remove_isolated_vertices(mesh)
    mesh, _ = pymesh.remove_duplicated_vertices(mesh, 0.001)

    return mesh


def compute_apbs(vertices, pdb_file, tmp_file_base):
    '''Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex'''

    args = ['pdb2pqr',
            '--ff=parse',
            '--whitespace',
            '--noopt',
            '--apbs-input',
            os.path.basename(pdb_file),
            os.path.basename(tmp_file_base)]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    args = ['apbs', os.path.basename(tmp_file_base + '.in')]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    with open(tmp_file_base + '.csv', 'w') as vertfile:
        for vert in vertices:
            vertfile.write(f'{vert[0]},{vert[1]},{vert[2]}\n')

    args = ['multivalue',
            os.path.basename(tmp_file_base) + '.csv',
            os.path.basename(tmp_file_base) + '.dx',
            os.path.basename(tmp_file_base) + '_out.csv']
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    with open(tmp_file_base + '_out.csv') as chargefile:
        charges = np.array([0.0] * len(vertices))
        for ix, line in enumerate(chargefile.readlines()):
            charges[ix] = float(line.split(',')[3])

    return charges


def compute_normal(vertex, face):
    '''Compute the normal of a triangulation

    Copyright (c) 2004 Gabriel Peyr
    Converted to Python by Pablo Gainza LPDI EPFL 2017

    Example:
        >>> normal,normalf = compute_normal(vertex, face)

    Args:
        vertex: 3xn matrix of vertices
        face: 3xm matrix of face indices.

    Returns:
      normal(i,:) is the normal at vertex i.
      normalf(j,:) is the normal at face j.

    '''

    vertex = vertex.T
    face = face.T

    nface = np.size(face, 1)
    nvert = np.size(vertex, 1)

    normal = np.zeros((3, nvert))

    # unit normals to the faces
    normalf = crossp(vertex[:, face[1, :]] - vertex[:, face[0, :]],
                     vertex[:, face[2, :]] - vertex[:, face[0, :]])

    sum_squares = np.sum(normalf ** 2, 0)

    d = np.sqrt(sum_squares)
    d[d < EPSILON] = 1

    normalf = normalf / repmat(d, 3, 1)

    # unit normal to the vertex
    normal = np.zeros((3, nvert))
    for i in np.arange(0, nface):
        f = face[:, i]

        for j in np.arange(3):
            normal[:, f[j]] = normal[:, f[j]] + normalf[:, i]

    # normalize
    d = np.sqrt(np.sum(normal ** 2, 0))
    d[d < EPSILON] = 1
    normal = normal / repmat(d, 3, 1)

    # enforce that the normal are outward
    vertex_means = np.mean(vertex, 0)

    v = vertex - repmat(vertex_means, 3, 1)
    s = np.sum(np.multiply(v, normal), 1)

    if np.sum(s > 0) < np.sum(s < 0):
        # flip
        normal = -normal
        normalf = -normalf

    return normal.T


def crossp(x, y):
    '''x and y are (m,3) dimensional'''
    z = np.zeros((x.shape))

    z[0, :] = np.multiply(x[1, :], y[2, :]) - np.multiply(x[2, :], y[1, :])
    z[1, :] = np.multiply(x[2, :], y[0, :]) - np.multiply(x[0, :], y[2, :])
    z[2, :] = np.multiply(x[0, :], y[1, :]) - np.multiply(x[1, :], y[0, :])

    return z
