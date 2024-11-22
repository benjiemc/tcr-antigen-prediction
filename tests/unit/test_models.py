import os
from unittest import TestCase

import numpy as np
import torch

from tcr_antigen_prediction.models import MasifPPISearch

DATA_DIR = 'tests/data'


class TestMasifPPISearch(TestCase):
    def test(self):
        model = MasifPPISearch(12.0, n_thetas=16, n_rhos=5, n_rotations=16)
        model.load_state_dict(torch.load(os.path.join(DATA_DIR, 'test_weights.pt'), weights_only=True))

        input_feats = np.load(os.path.join(DATA_DIR, 'input_feats.npy')).astype(np.float32)
        rho_coords = np.load(os.path.join(DATA_DIR, 'rho_wrt_center.npy')).astype(np.float32)
        theta_coords = np.load(os.path.join(DATA_DIR, 'theta_wrt_center.npy')).astype(np.float32)
        mask = np.load(os.path.join(DATA_DIR, 'mask.npy')).astype(np.float32)

        rho_coords = np.expand_dims(rho_coords, 2)
        theta_coords = np.expand_dims(theta_coords, 2)
        mask = np.expand_dims(mask, 2)

        input_feats = torch.from_numpy(input_feats)
        rho_coords = torch.from_numpy(rho_coords)
        theta_coords = torch.from_numpy(theta_coords)
        mask = torch.from_numpy(mask)

        import pdb; pdb.set_trace()
        output = model(input_feats, rho_coords, theta_coords, mask)
        np.testing.assert_array_almost_equal(output.detach().numpy(), np.array([]))
