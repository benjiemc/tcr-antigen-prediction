import os
from unittest import TestCase

import numpy as np
from Bio.PDB import PDBParser

from tcr_antigen_prediction.comparisons import compute_structural_distances, rmsd

TEST_DATA = 'tests/data'


class TestRmsd(TestCase):
    def test(self):
        arr1 = np.array([[1.00, 0.00, 0.00],
                         [3.00, 0.00, 1.00],
                         [2.00, 0.00, 1.00]])
        arr2 = np.array([[6.00, 0.00, 0.00],
                         [1.00, 0.70, 1.00],
                         [4.00, 0.00, 1.00]])

        self.assertAlmostEqual(rmsd(arr1, arr2), 3.3411575)

    def test_zeros(self):
        arr1 = np.zeros((10, 3))
        arr2 = np.zeros((10, 3))

        self.assertAlmostEqual(rmsd(arr1, arr2), 0.00)


class ComputeStructuralDistances(TestCase):
    def test(self):
        pdb_parser = PDBParser(QUIET=True)
        struct_decab = pdb_parser.get_structure('decab', os.path.join(TEST_DATA, '7q9b_DECAB.pdb'))
        struct_ijhfg = pdb_parser.get_structure('ijhfg', os.path.join(TEST_DATA, '7q9b_IJHFG.pdb'))

        distance_matrix = compute_structural_distances([struct_decab, struct_ijhfg],
                                                       [{'alpha_chain': 'D',
                                                         'beta_chain': 'E',
                                                         'antigen_chain': 'C',
                                                         'mhc_chain1': 'A',
                                                         'mhc_chain2': 'B'},
                                                        {'alpha_chain': 'I',
                                                         'beta_chain': 'J',
                                                         'antigen_chain': 'H',
                                                         'mhc_chain1': 'F',
                                                         'mhc_chain2': 'G'}],
                                                       'MH1')

        np.testing.assert_array_almost_equal(distance_matrix, [[0.0, 1.152162],
                                                               [1.152162, 0.0]])
