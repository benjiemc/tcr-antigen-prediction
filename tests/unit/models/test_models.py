from unittest import TestCase

import numpy as np
import torch

from tcr_antigen_prediction.models import TCRContactMapPredictor


class TestTCRContactMapPredictor(TestCase):
    def test(self):
        torch.manual_seed(0)

        cdr_1_peptide_contact_map = torch.tensor(
            [
                [0.03, 0.03, 0.03, 0.03, 0.03],
                [0.05, 0.09, 0.09, 0.09, 0.03],
                [0.05, 0.09, 0.09, 0.09, 0.03],
                [0.03, 0.03, 0.03, 0.03, 0.03],
            ]
        )
        cdr_2_peptide_contact_map = torch.tensor(
            [
                [0.03, 0.03, 0.03, 0.03, 0.03],
                [0.05, 0.09, 0.09, 0.09, 0.03],
                [0.05, 0.09, 0.09, 0.09, 0.03],
                [0.03, 0.03, 0.03, 0.03, 0.03],
            ]
        )
        cdr_3_peptide_contact_map = torch.tensor(
            [
                [0.03, 0.03, 0.03, 0.03, 0.03],
                [0.03, 0.05, 0.05, 0.03, 0.03],
                [0.03, 0.05, 0.05, 0.03, 0.03],
                [0.03, 0.03, 0.05, 0.03, 0.03],
                [0.03, 0.03, 0.03, 0.03, 0.03],
                [0.03, 0.03, 0.03, 0.03, 0.03],
            ]
        )

        cdr_peptide_contact_maps = (
            cdr_1_peptide_contact_map,
            cdr_2_peptide_contact_map,
            cdr_3_peptide_contact_map,
            cdr_1_peptide_contact_map,
            cdr_2_peptide_contact_map,
            cdr_3_peptide_contact_map,
        )

        cdr_1_mhc_contact_map = torch.tensor(
            [
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.05, 0.09, 0.09, 0.09, 0.03, 0.0, 0.0],
                [0.0, 0.05, 0.09, 0.09, 0.09, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
            ]
        )
        cdr_2_mhc_contact_map = torch.tensor(
            [
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.05, 0.09, 0.09, 0.09, 0.03, 0.0, 0.0],
                [0.0, 0.05, 0.09, 0.09, 0.09, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
            ]
        )
        cdr_3_mhc_contact_map = torch.tensor(
            [
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.05, 0.05, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.05, 0.05, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.03, 0.05, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
                [0.0, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0, 0.0],
            ]
        )

        cdr_mhc_contact_maps = (
            cdr_1_mhc_contact_map,
            cdr_2_mhc_contact_map,
            cdr_3_mhc_contact_map,
            cdr_1_mhc_contact_map,
            cdr_2_mhc_contact_map,
            cdr_3_mhc_contact_map,
        )

        model = TCRContactMapPredictor(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            cdr_1_length=4,
            cdr_2_length=4,
            cdr_3_length=6,
            peptide_length=5,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr_1a_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr_2a_batch = cdr_1a_batch.clone().detach()
        cdr_1b_batch = cdr_1a_batch.clone().detach()
        cdr_2b_batch = cdr_1a_batch.clone().detach()

        cdr_3a_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr_3b_batch = cdr_3a_batch.clone().detach()

        peptide_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
            ],
            dtype=torch.float32,
        )

        mhc_batch = torch.tensor(
            [
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                ],
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                ],
            ],
            dtype=torch.float32,
        )

        predictions = model.forward(
            cdr_1a_batch,
            cdr_2a_batch,
            cdr_3a_batch,
            cdr_1b_batch,
            cdr_2b_batch,
            cdr_3b_batch,
            peptide_batch,
            mhc_batch,
        )
        np.testing.assert_array_almost_equal(
            predictions.detach(), np.array([[0.027157, -0.008714], [0.02683, -0.008769]])
        )

    def test_no_contact_maps(self):
        torch.manual_seed(0)

        model = TCRContactMapPredictor(
            cdr_1_length=4,
            cdr_2_length=4,
            cdr_3_length=6,
            peptide_length=5,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr_1a_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr_2a_batch = cdr_1a_batch.clone().detach()
        cdr_1b_batch = cdr_1a_batch.clone().detach()
        cdr_2b_batch = cdr_1a_batch.clone().detach()

        cdr_3a_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr_3b_batch = cdr_3a_batch.clone().detach()

        peptide_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
            ],
            dtype=torch.float32,
        )

        mhc_batch = torch.tensor(
            [
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                ],
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0],
                ],
            ],
            dtype=torch.float32,
        )

        predictions = model.forward(
            cdr_1a_batch,
            cdr_2a_batch,
            cdr_3a_batch,
            cdr_1b_batch,
            cdr_2b_batch,
            cdr_3b_batch,
            peptide_batch,
            mhc_batch,
        )
        np.testing.assert_array_almost_equal(
            predictions.detach(), np.array([[0.026841, -0.008762], [0.026573, -0.008937]])
        )
