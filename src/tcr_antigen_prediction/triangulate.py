import os
from subprocess import Popen, PIPE

import numpy as np
import pymesh
from Bio.PDB import PDBParser, Selection, NeighborSearch
from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral
from numpy.matlib import repmat
from sklearn.neighbors import KDTree

from tcr_antigen_prediction.io import output_pdb_as_xyzrn, read_msms
from tcr_antigen_prediction.chemistry import (donorAtom,
                                              polarHydrogens,
                                              acceptorAngleAtom,
                                              acceptorPlaneAtom,
                                              hbond_std_dev)

EPSILON = 1.0e-6


def computeCharges(pdb_filename, vertices, names):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_filename, pdb_filename)
    residues = {}
    for res in struct.get_residues():
        chain_id = res.get_parent().get_id()
        if chain_id == '':
            chain_id = ' '
        residues[(chain_id, res.get_id())] = res

    atoms = Selection.unfold_entities(struct, 'A')
    satisfied_CO, satisfied_HN = computeSatisfied_CO_HN(atoms)

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
        if atom_name == 'H' and res_id in satisfied_HN:
            continue
        if atom_name == 'O' and res_id in satisfied_CO:
            continue
        # Compute the charge of the vertex
        charge[ix] = computeChargeHelper(
            atom_name, residues[(chain_id, res_id)], vertices[ix]
        )

    return charge


# Compute the charge of a vertex in a residue.
def computeChargeHelper(atom_name, res, v):
    # Check if it is a polar hydrogen.
    if isPolarHydrogen(atom_name, res):
        donor_atom_name = donorAtom[atom_name]
        a = res[donor_atom_name].get_coord()  # N/O
        b = res[atom_name].get_coord()  # H
        # Donor-H is always 180.0 degrees, = pi
        angle_deviation = computeAngleDeviation(a, b, v, np.pi)
        angle_penalty = computeAnglePenalty(angle_deviation)
        return 1.0 * angle_penalty
    # Check if it is an acceptor oxygen or nitrogen
    elif isAcceptorAtom(atom_name, res):
        acceptor_atom = res[atom_name]
        b = acceptor_atom.get_coord()
        # try:
        a = res[acceptorAngleAtom[atom_name]].get_coord()
        # except:
        # return 0.0
        # 120 degress for acceptor
        angle_deviation = computeAngleDeviation(a, b, v, 2 * np.pi / 3)
        # TODO: This should not be 120 for all atoms, i.e. for HIS it should be
        #       ~125.0
        angle_penalty = computeAnglePenalty(angle_deviation)
        plane_penalty = 1.0
        if atom_name in acceptorPlaneAtom:
            # try:
            d = res[acceptorPlaneAtom[atom_name]].get_coord()
            # except:
            # return 0.0
            plane_deviation = computePlaneDeviation(d, a, b, v)
            plane_penalty = computeAnglePenalty(plane_deviation)
        return -1.0 * angle_penalty * plane_penalty
        # Compute the
    return 0.0


# Compute the absolute value of the deviation from theta
def computeAngleDeviation(a, b, c, theta):
    return abs(calc_angle(Vector(a), Vector(b), Vector(c)) - theta)


# Compute the angle deviation from a plane
def computePlaneDeviation(a, b, c, d):
    dih = calc_dihedral(Vector(a), Vector(b), Vector(c), Vector(d))
    dev1 = abs(dih)
    dev2 = np.pi - abs(dih)
    return min(dev1, dev2)


# angle_deviation from ideal value. TODO: do a more data-based solution
def computeAnglePenalty(angle_deviation):
    # Standard deviation: hbond_std_dev
    return max(0.0, 1.0 - (angle_deviation / (hbond_std_dev)) ** 2)


def isPolarHydrogen(atom_name, res):
    if atom_name in polarHydrogens[res.get_resname()]:
        return True
    else:
        return False


def isAcceptorAtom(atom_name, res):
    if atom_name.startswith('O'):
        return True
    else:
        if res.get_resname() == 'HIS':
            if atom_name == 'ND1' and 'HD1' not in res:
                return True
            if atom_name == 'NE2' and 'HE2' not in res:
                return True
    return False


# Compute the list of backbone C=O:H-N that are satisfied. These will be ignored.
def computeSatisfied_CO_HN(atoms):
    ns = NeighborSearch(atoms)
    satisfied_CO = set()
    satisfied_HN = set()
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
                        angle_N_H_O_dev = computeAngleDeviation(
                            res2['N'].get_coord(),
                            atom2.get_coord(),
                            atom1.get_coord(),
                            np.pi,
                        )
                        # Compute angle H:O=C, ideal value is ~160 +- 20 = 8*pi/9
                        angle_H_O_C_dev = computeAngleDeviation(
                            atom2.get_coord(),
                            atom1.get_coord(),
                            res1['C'].get_coord(),
                            8 * np.pi / 9,
                        )
                        # Allowed deviations: 30 degrees (pi/6) and 20 degrees
                        #       (pi/9)
                        if angle_N_H_O_dev - np.pi / 6 < 0 and angle_H_O_C_dev - np.pi / 9 < 0.0:
                            satisfied_CO.add(res1.get_id())
                            satisfied_HN.add(res2.get_id())
    return satisfied_CO, satisfied_HN


def assignChargesToNewMesh(new_vertices, old_vertices, old_charges):
    '''
    Compute the charge of a new mesh, based on the charge of an old mesh.
    Use the top vertex in distance, for now (later this should be smoothed over 3
    or 4 vertices)
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
            new_charges[vi_new] += (
                old_charges[vi_old[i]] * (1 / dist_old[i]) / total_dist
            )
    return new_charges


def computeMSMS(pdb_file,  protonate=True):
    file_base = pdb_file.rsplit('.', 1)[0]
    out_xyzrn = file_base + '.xyzrn'

    if protonate:
        output_pdb_as_xyzrn(pdb_file, out_xyzrn)

    # Now run MSMS on xyzrn file
    args = ['msms',
            '-density',
            '3.0',
            '-hdensity',
            '3.0',
            '-probe',
            '1.5',
            '-if', out_xyzrn,
            '-of', file_base,
            '-af', file_base]

    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    _, _ = p2.communicate()

    vertices, faces, normals, names = read_msms(file_base)
    areas = {}
    ses_file = open(file_base + '.area')
    next(ses_file)  # ignore header line

    for line in ses_file:
        fields = line.split()
        areas[fields[3]] = fields[1]

    # Remove temporary files.
    # os.remove(file_base+'.area')
    # os.remove(file_base+'.xyzrn')
    # os.remove(file_base+'.vert')
    # os.remove(file_base+'.face')

    return vertices, faces, normals, names, areas


# Kyte Doolittle scale
kd_scale = {}
kd_scale['ILE'] = 4.5
kd_scale['VAL'] = 4.2
kd_scale['LEU'] = 3.8
kd_scale['PHE'] = 2.8
kd_scale['CYS'] = 2.5
kd_scale['MET'] = 1.9
kd_scale['ALA'] = 1.8
kd_scale['GLY'] = -0.4
kd_scale['THR'] = -0.7
kd_scale['SER'] = -0.8
kd_scale['TRP'] = -0.9
kd_scale['TYR'] = -1.3
kd_scale['PRO'] = -1.6
kd_scale['HIS'] = -3.2
kd_scale['GLU'] = -3.5
kd_scale['GLN'] = -3.5
kd_scale['ASP'] = -3.5
kd_scale['ASN'] = -3.5
kd_scale['LYS'] = -3.9
kd_scale['ARG'] = -4.5


def computeHydrophobicity(names):
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
    print('Removing degenerated triangles')
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


def computeAPBS(vertices, pdb_file, tmp_file_base):
    '''
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    '''
    args = [
        'pdb2pqr',
        '--ff=parse',
        '--whitespace',
        '--noopt',
        '--apbs-input',
        os.path.basename(pdb_file),
        os.path.basename(tmp_file_base),
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    args = ['apbs', os.path.basename(tmp_file_base + '.in')]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    vertfile = open(tmp_file_base + '.csv', 'w')
    for vert in vertices:
        vertfile.write('{},{},{}\n'.format(vert[0], vert[1], vert[2]))
    vertfile.close()

    args = [
        'multivalue',
        os.path.basename(tmp_file_base) + '.csv',
        os.path.basename(tmp_file_base) + '.dx',
        os.path.basename(tmp_file_base) + '_out.csv',
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(tmp_file_base))
    _, _ = p2.communicate()

    # Read the charge file
    chargefile = open(tmp_file_base + '_out.csv')
    charges = np.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        charges[ix] = float(line.split(',')[3])

    return charges


def compute_normal(vertex, face):
    '''
    compute_normal - compute the normal of a triangulation
    vertex: 3xn matrix of vertices
    face: 3xm matrix of face indices.

      normal,normalf = compute_normal(vertex,face)

      normal(i,:) is the normal at vertex i.
      normalf(j,:) is the normal at face j.

    Copyright (c) 2004 Gabriel Peyr
    Converted to Python by Pablo Gainza LPDI EPFL 2017
    '''

    vertex = vertex.T
    face = face.T
    nface = np.size(face, 1)
    nvert = np.size(vertex, 1)
    normal = np.zeros((3, nvert))
    # unit normals to the faces
    normalf = crossp(
        vertex[:, face[1, :]] - vertex[:, face[0, :]],
        vertex[:, face[2, :]] - vertex[:, face[0, :]],
    )
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
    ''' x and y are (m,3) dimensional '''
    z = np.zeros((x.shape))
    z[0, :] = np.multiply(x[1, :], y[2, :]) - np.multiply(x[2, :], y[1, :])
    z[1, :] = np.multiply(x[2, :], y[0, :]) - np.multiply(x[0, :], y[2, :])
    z[2, :] = np.multiply(x[0, :], y[1, :]) - np.multiply(x[1, :], y[0, :])
    return z
