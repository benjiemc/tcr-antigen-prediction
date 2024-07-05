'''Functions for inputting and outputing various file formats.'''
import pymesh
import numpy as np
from Bio.PDB import PDBParser

from tcr_antigen_prediction.chemistry import RADII, POLAR_HYDROGENS
from tcr_antigen_prediction.geometry import compute_polar_coordinates
from tcr_antigen_prediction.surface import compute_ddc, normalize_electrostatics


def output_pdb_as_xyzrn(pdb_filename, xyzrn_filename):
    '''Write PDB file as xyzrn.

    Args:
        pdb_filename: input pdb filename
        xyzrn_filename: output in xyzrn format.

    '''
    parser = PDBParser()
    struct = parser.get_structure(pdb_filename, pdb_filename)

    with open(xyzrn_filename, 'w') as outfile:
        for atom in struct.get_atoms():
            name = atom.get_name()
            residue = atom.get_parent()

            # Ignore hetatms.
            if residue.get_id()[0] != ' ':
                continue

            resname = residue.get_resname()
            chain = residue.get_parent().get_id()
            atomtype = name[0]

            color = 'Green'
            coords = None

            if atomtype in RADII and resname in POLAR_HYDROGENS:
                if atomtype == 'O':
                    color = 'Red'

                if atomtype == 'N':
                    color = 'Blue'

                if atomtype == 'H':
                    if name in POLAR_HYDROGENS[resname]:
                        color = 'Blue'  # Polar hydrogens

                coords = f'{atom.get_coord()[0]: .06f} {atom.get_coord()[1]: .06f} {atom.get_coord()[2]: .06f}'
                insertion = 'x'

                if residue.get_id()[2] != ' ':
                    insertion = residue.get_id()[2]

                full_id = f'{chain}_{residue.get_id()[1]}_{insertion}_{resname}_{name}_{color}'

            if coords is not None:
                outfile.write(coords + ' ' + RADII[atomtype] + ' 1 ' + full_id + '\n')


def read_data_from_surface(ply_fn, max_distance, max_shape_size):
    '''Read data from a ply file and decompose into patches.

    Returns:
        A tuple containing:
          - list of features per patch
          - list of angular and polar coordinates.
          - list of indices of neighbors in the patch.
          - list of shape complementarity labels (computed here).

    '''
    mesh = pymesh.load_mesh(ply_fn)

    # Normals:
    n1 = mesh.get_attribute('vertex_nx')
    n2 = mesh.get_attribute('vertex_ny')
    n3 = mesh.get_attribute('vertex_nz')

    normals = np.stack([n1, n2, n3], axis=1)

    # Compute the angular and radial coordinates.
    rho, theta, neigh_indices, mask = compute_polar_coordinates(mesh,
                                                                radius=max_distance,
                                                                max_vertices=max_shape_size)

    # Compute the principal curvature components for the shape index.
    mesh.add_attribute('vertex_mean_curvature')
    h_mat = mesh.get_attribute('vertex_mean_curvature')
    mesh.add_attribute('vertex_gaussian_curvature')
    k_mat = mesh.get_attribute('vertex_gaussian_curvature')

    elem = np.square(h_mat) - k_mat

    # In some cases this equation is less than zero, likely due to the method that computes the mean and gaussian
    # curvature. set to an epsilon.
    elem[elem < 0] = 1e-8

    k1 = h_mat + np.sqrt(elem)
    k2 = h_mat - np.sqrt(elem)

    # Compute the shape index
    si = (k1 + k2) / (k1 - k2)
    si = np.arctan(si) * (2 / np.pi)

    # Normalize the charge.
    charge = mesh.get_attribute('vertex_charge')
    charge = normalize_electrostatics(charge)

    # Hbond features
    hbond = mesh.get_attribute('vertex_hbond')

    # Hydropathy features
    # Normalize hydropathy by dividing by 4.5
    hphob = mesh.get_attribute('vertex_hphob') / 4.5

    # Iface labels (for ground truth only)
    if 'vertex_iface' in mesh.get_attribute_names():
        iface_labels = mesh.get_attribute('vertex_iface')

    else:
        iface_labels = np.zeros_like(hphob)

    # n: number of patches, equal to the number of vertices.
    n = len(mesh.vertices)

    input_feat = np.zeros((n, max_shape_size, 5))

    # Compute the input features for each patch.
    for vix in range(n):
        # Patch members.
        neigh_vix = np.array(neigh_indices[vix])

        # Compute the distance-dependent curvature for all neighbors of the patch.
        patch_v = mesh.vertices[neigh_vix]
        patch_n = normals[neigh_vix]
        patch_cp = np.where(neigh_vix == vix)[0][0]  # central point

        mask_pos = np.where(mask[vix] == 1.0)[0]  # nonzero elements
        patch_rho = rho[vix][mask_pos]  # nonzero elements of rho

        ddc = compute_ddc(patch_v, patch_n, patch_cp, patch_rho)

        input_feat[vix, :len(neigh_vix), 0] = si[neigh_vix]
        input_feat[vix, :len(neigh_vix), 1] = ddc
        input_feat[vix, :len(neigh_vix), 2] = hbond[neigh_vix]
        input_feat[vix, :len(neigh_vix), 3] = charge[neigh_vix]
        input_feat[vix, :len(neigh_vix), 4] = hphob[neigh_vix]

    return input_feat, rho, theta, mask, neigh_indices, iface_labels, np.copy(mesh.vertices)


def read_msms(file_root):
    '''read the surface from the msms output. MSMS outputs two files: {file_root}.vert and {file_root}.face'''
    with open(file_root + '.vert') as vertfile:
        meshdata = vertfile.read().rstrip().split('\n')

    # Read number of vertices.
    count = {}
    header = meshdata[2].split()
    count['vertices'] = int(header[0])

    # Data Structures
    vertices = np.zeros((count['vertices'], 3))
    normalv = np.zeros((count['vertices'], 3))

    atom_id = [''] * count['vertices']
    res_id = [''] * count['vertices']

    for i in range(3, len(meshdata)):
        fields = meshdata[i].split()
        vi = i - 3

        vertices[vi][0] = float(fields[0])
        vertices[vi][1] = float(fields[1])
        vertices[vi][2] = float(fields[2])

        normalv[vi][0] = float(fields[3])
        normalv[vi][1] = float(fields[4])
        normalv[vi][2] = float(fields[5])

        atom_id[vi] = fields[7]
        res_id[vi] = fields[9]

        count['vertices'] -= 1

    with open(file_root + '.face') as facefile:
        meshdata = facefile.read().rstrip().split('\n')

    # Read number of vertices.
    header = meshdata[2].split()
    count['faces'] = int(header[0])
    faces = np.zeros((count['faces'], 3), dtype=int)

    for i in range(3, len(meshdata)):
        fi = i - 3
        fields = meshdata[i].split()

        faces[fi][0] = int(fields[0]) - 1
        faces[fi][1] = int(fields[1]) - 1
        faces[fi][2] = int(fields[2]) - 1

        count['faces'] -= 1

    assert count['vertices'] == 0
    assert count['faces'] == 0

    return vertices, faces, normalv, res_id


def save_ply(filename,
             vertices,
             faces=None,
             normals=None,
             charges=None,
             vertex_cb=None,
             hbond=None,
             hphob=None,
             iface=None,
             normalize_charges=False):
    ''' Save vertices, mesh in ply format.

    Args:
        vertices: coordinates of vertices
        faces: mesh

    '''
    if faces is None:
        faces = []

    mesh = pymesh.form_mesh(vertices, faces)

    if normals is not None:
        n1 = normals[:, 0]
        n2 = normals[:, 1]
        n3 = normals[:, 2]

        mesh.add_attribute('vertex_nx')
        mesh.set_attribute('vertex_nx', n1)

        mesh.add_attribute('vertex_ny')
        mesh.set_attribute('vertex_ny', n2)

        mesh.add_attribute('vertex_nz')
        mesh.set_attribute('vertex_nz', n3)

    if charges is not None:
        mesh.add_attribute('charge')

        if normalize_charges:
            charges = charges / 10

        mesh.set_attribute('charge', charges)

    if hbond is not None:
        mesh.add_attribute('hbond')
        mesh.set_attribute('hbond', hbond)

    if vertex_cb is not None:
        mesh.add_attribute('vertex_cb')
        mesh.set_attribute('vertex_cb', vertex_cb)

    if hphob is not None:
        mesh.add_attribute('vertex_hphob')
        mesh.set_attribute('vertex_hphob', hphob)

    if iface is not None:
        mesh.add_attribute('vertex_iface')
        mesh.set_attribute('vertex_iface', iface)

    pymesh.save_mesh(filename, mesh, *mesh.get_attribute_names(), use_float=True, ascii=True)
