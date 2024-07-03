import argparse
import logging
import time
import os
import numpy as np

from tcr_antigen_prediction.apps._log import setup_logger
from tcr_antigen_prediction.read_data_from_surface import read_data_from_surface, compute_shape_complementarity

logger = logging.getLogger()

np.random.seed(0)

parser = argparse.ArgumentParser()

parser.add_argument('--mode', choices=['ppi_search', 'site'], required=True,
                    help='mode to compute features')
parser.add_argument('--output', '-o', required=True, help='output path')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")
parser.add_argument('input', nargs=2, help='path to input ply files')

config = {}

config["ppi_search"] = {}
config["ppi_search"]["max_shape_size"] = 200
config["ppi_search"]["max_distance"] = 12.0
# Parameters for shape complementarity calculations.
config["ppi_search"]["sc_radius"] = 12.0
config["ppi_search"]["sc_interaction_cutoff"] = 1.5
config["ppi_search"]["sc_w"] = 0.25

config["site"] = {}
config["site"]["max_shape_size"] = 100
config["site"]["max_distance"] = 9.0


if __name__ == '__main__':
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    total_shapes = 0
    total_ppi_pairs = 0
    np.random.seed(0)

    all_list_desc = []
    all_list_coords = []
    all_list_shape_idx = []
    all_list_names = []
    idx_positives = []
    print('Reading data from input ply surface files.')

    pids = [f'p{num}' for num in range(1, len(args.input) + 1)]

    # Compute shape complementarity between the two proteins.
    rho = {}
    neigh_indices = {}
    mask = {}
    input_feat = {}
    theta = {}
    iface_labels = {}
    verts = {}

    for pid, ply_file in zip(pids, args.input):
        (input_feat[pid],
         rho[pid],
         theta[pid],
         mask[pid],
         neigh_indices[pid],
         iface_labels[pid],
         verts[pid]) = read_data_from_surface(ply_file,
                                              max_distance=config[args.mode]['max_distance'],
                                              max_shape_size=config[args.mode]['max_shape_size'])

    if len(pids) > 1 and args.mode == 'masif_ppi_search':
        start_time = time.time()
        p1_sc_labels, p2_sc_labels = compute_shape_complementarity(
            ply_file['p1'], ply_file['p2'],
            neigh_indices['p1'], neigh_indices['p2'],
            rho['p1'], rho['p2'],
            mask['p1'], mask['p2'],
            sc_w=config[args.mode]['sc_w'],
            sc_interaction_cutoff=config[args.mode]['sc_interaction_count'],
            sc_radius=config[args.mode]['sc_radius'],
        )
        np.save(os.path.join(args.output, 'p1_sc_labels'), p1_sc_labels)
        np.save(os.path.join(args.output, 'p2_sc_labels'), p2_sc_labels)
        end_time = time.time()
        print("Computing shape complementarity took {:.2f}".format(end_time - start_time))

    # Save data only if everything went well.
    for pid in pids:
        np.save(os.path.join(args.output, pid + '_rho_wrt_center'), rho[pid])
        np.save(os.path.join(args.output, pid + '_theta_wrt_center'), theta[pid])
        np.save(os.path.join(args.output, pid + '_input_feat'), input_feat[pid])
        np.save(os.path.join(args.output, pid + '_mask'), mask[pid])
        np.save(os.path.join(args.output, pid + '_list_indices'), neigh_indices[pid])
        np.save(os.path.join(args.output, pid + '_iface_labels'), iface_labels[pid])

        # Save x, y, z
        np.save(os.path.join(args.output, pid + '_X.npy'), verts[pid][:, 0])
        np.save(os.path.join(args.output, pid + '_Y.npy'), verts[pid][:, 1])
        np.save(os.path.join(args.output, pid + '_Z.npy'), verts[pid][:, 2])
