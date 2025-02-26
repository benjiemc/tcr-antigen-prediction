from unittest import TestCase

import numpy as np
import torch

from tcr_antigen_prediction.models import TCRStructMap


class TestTCRStructMap(TestCase):
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

        model = TCRStructMap(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            cdr1_alpha_length=4,
            cdr2_alpha_length=4,
            cdr3_alpha_length=6,
            cdr1_beta_length=4,
            cdr2_beta_length=4,
            cdr3_beta_length=6,
            peptide_length=5,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr1_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr2_alpha_batch = cdr1_alpha_batch.clone().detach()
        cdr1_beta_batch = cdr1_alpha_batch.clone().detach()
        cdr2_beta_batch = cdr1_alpha_batch.clone().detach()

        cdr3_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr3_beta_batch = cdr3_alpha_batch.clone().detach()

        peptide_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
            ],
            dtype=torch.float32,
        )

        mhc_pseudo_batch = torch.tensor(
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
            cdr1_alpha_batch,
            cdr2_alpha_batch,
            cdr3_alpha_batch,
            cdr1_beta_batch,
            cdr2_beta_batch,
            cdr3_beta_batch,
            peptide_batch,
            mhc_pseudo_batch,
        )
        np.testing.assert_array_almost_equal(predictions.detach(), np.array([[0.491795], [0.491714]]))

    def test_no_contact_maps(self):
        torch.manual_seed(0)

        model = TCRStructMap(
            cdr1_alpha_length=4,
            cdr2_alpha_length=4,
            cdr3_alpha_length=6,
            cdr1_beta_length=4,
            cdr2_beta_length=4,
            cdr3_beta_length=6,
            peptide_length=5,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr1_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr2_alpha_batch = cdr1_alpha_batch.clone().detach()
        cdr1_beta_batch = cdr1_alpha_batch.clone().detach()
        cdr2_beta_batch = cdr1_alpha_batch.clone().detach()

        cdr3_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr3_beta_batch = cdr3_alpha_batch.clone().detach()

        peptide_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
            ],
            dtype=torch.float32,
        )

        mhc_pseudo_batch = torch.tensor(
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
            cdr1_alpha_batch,
            cdr2_alpha_batch,
            cdr3_alpha_batch,
            cdr1_beta_batch,
            cdr2_beta_batch,
            cdr3_beta_batch,
            peptide_batch,
            mhc_pseudo_batch,
        )
        np.testing.assert_array_almost_equal(predictions.detach(), np.array([[0.491717], [0.491649]]))

    def test_no_peptide(self):
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

        model = TCRStructMap(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            cdr1_alpha_length=4,
            cdr2_alpha_length=4,
            cdr3_alpha_length=6,
            cdr1_beta_length=4,
            cdr2_beta_length=4,
            cdr3_beta_length=6,
            peptide_length=0,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr1_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr2_alpha_batch = cdr1_alpha_batch.clone().detach()
        cdr1_beta_batch = cdr1_alpha_batch.clone().detach()
        cdr2_beta_batch = cdr1_alpha_batch.clone().detach()

        cdr3_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr3_beta_batch = cdr3_alpha_batch.clone().detach()

        mhc_pseudo_batch = torch.tensor(
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
            cdr1_alpha_batch,
            cdr2_alpha_batch,
            cdr3_alpha_batch,
            cdr1_beta_batch,
            cdr2_beta_batch,
            cdr3_beta_batch,
            None,
            mhc_pseudo_batch,
        )
        np.testing.assert_array_almost_equal(predictions.detach(), np.array([[0.510418], [0.51038]]))

    def test_no_peptide_no_contact_maps(self):
        torch.manual_seed(0)

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

        model = TCRStructMap(
            None,
            cdr_mhc_contact_maps,
            cdr1_alpha_length=4,
            cdr2_alpha_length=4,
            cdr3_alpha_length=6,
            cdr1_beta_length=4,
            cdr2_beta_length=4,
            cdr3_beta_length=6,
            peptide_length=0,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr1_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr2_alpha_batch = cdr1_alpha_batch.clone().detach()
        cdr1_beta_batch = cdr1_alpha_batch.clone().detach()
        cdr2_beta_batch = cdr1_alpha_batch.clone().detach()

        cdr3_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr3_beta_batch = cdr3_alpha_batch.clone().detach()

        mhc_pseudo_batch = torch.tensor(
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
            cdr1_alpha_batch,
            cdr2_alpha_batch,
            cdr3_alpha_batch,
            cdr1_beta_batch,
            cdr2_beta_batch,
            cdr3_beta_batch,
            None,
            mhc_pseudo_batch,
        )
        np.testing.assert_array_almost_equal(predictions.detach(), np.array([[0.510418], [0.51038]]))

    def test_no_cdr_1_or_2(self):
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

        model = TCRStructMap(
            cdr_peptide_contact_maps,
            cdr_mhc_contact_maps,
            cdr1_alpha_length=0,
            cdr2_alpha_length=0,
            cdr3_alpha_length=6,
            cdr1_beta_length=0,
            cdr2_beta_length=0,
            cdr3_beta_length=6,
            peptide_length=5,
            mhc_length=8,
            drop_out_rate=0.0,
            input_depth=4,
        )

        cdr3_alpha_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
            ],
            dtype=torch.float32,
        )

        cdr3_beta_batch = cdr3_alpha_batch.clone().detach()

        peptide_batch = torch.tensor(
            [
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
            ],
            dtype=torch.float32,
        )

        mhc_pseudo_batch = torch.tensor(
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
            None,
            None,
            cdr3_alpha_batch,
            None,
            None,
            cdr3_beta_batch,
            peptide_batch,
            mhc_pseudo_batch,
        )
        np.testing.assert_array_almost_equal(predictions.detach(), np.array([[0.518837], [0.518778]]))
