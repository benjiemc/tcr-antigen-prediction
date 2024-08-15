from unittest import TestCase
import os

import numpy as np

from tcr_antigen_prediction.triangulate import compute_msms

TEST_DATA = 'tests/data'


class TestComputeMSMS(TestCase):
    def test(self):
        vertices, faces, normals, names, areas = compute_msms(os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

        np.testing.assert_array_almost_equal(vertices[:3], np.array([[49.032, -8.53, 46.926],
                                                                    [49.531, -8.33, 47.252],
                                                                    [50.124, -8.236, 47.061]]))
        np.testing.assert_array_almost_equal(vertices[-3:], np.array([[80.817, -16.696, 17.894],
                                                                    [80.817, -16.831, 17.302],
                                                                    [80.817, -16.696, 16.71]]))

        np.testing.assert_array_equal(faces[:3], np.array([[0, 12238, 1],
                                                        [1, 12238, 12239],
                                                        [12238, 12240, 12239]]))
        np.testing.assert_array_equal(faces[-3:], np.array([[12175, 12190, 12174],
                                                            [12187, 12183, 12184],
                                                            [12178, 12181, 12177]]))

        np.testing.assert_array_almost_equal(normals[:3], np.array([[0.222,  0.858, -0.463],
                                                                    [-0.11 ,  0.725, -0.68 ],
                                                                    [-0.506,  0.662, -0.552]]))
        np.testing.assert_array_almost_equal(normals[-3:], np.array([[-0.223, -0.878, 0.423],
                                                                    [-0.223, -0.975, 0.00],
                                                                    [-0.223, -0.878, -0.423]]))

        np.testing.assert_array_equal(names[:3], np.array(['A_17_x_ARG_NH2_Blue',
                                                        'A_17_x_ARG_CZ_Green',
                                                        'A_17_x_ARG_NH1_Blue']))
        np.testing.assert_array_equal(names[-3:], np.array(['A_1087_x_GLU_OE2_Red',
                                                            'A_1087_x_GLU_OE2_Red',
                                                            'A_1087_x_GLU_OE2_Red']))

        self.assertEqual(areas['A_17_x_ARG_NH2_Blue'], '19.4041')
        self.assertEqual(areas['A_17_x_ARG_CZ_Green'], '8.5787')
        self.assertEqual(areas['A_17_x_ARG_NH1_Blue'], '20.8036')
        self.assertEqual(areas['A_1087_x_GLU_OE2_Red'], '14.7990')
