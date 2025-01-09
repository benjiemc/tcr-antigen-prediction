from unittest import TestCase

from tcr_antigen_prediction.data.imgt_numbering import assign_helix


class TestAssignHelix(TestCase):
    def test_mhc_i_alpha_helix(self):
        self.assertEqual(assign_helix('MH1', 'mhc_chain1', '30'), 'alpha')

    def test_mhc_i_beta_helix(self):
        self.assertEqual(assign_helix('MH1', 'mhc_chain1', '1030'), 'beta')

    def test_mhc_ii_alpha_helix(self):
        self.assertEqual(assign_helix('MH2', 'mhc_chain1', '30'), 'alpha')

    def test_mhc_ii_beta_helix(self):
        self.assertEqual(assign_helix('MH2', 'mhc_chain2', '30'), 'beta')
