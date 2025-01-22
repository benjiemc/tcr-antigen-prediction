from unittest import TestCase

from tcr_antigen_prediction.data.utils import centre_pad, mhc_code_to_slug, mhc_slug_to_code


class TestMHCCodeToSlug(TestCase):
    def test_human(self):
        self.assertEqual(mhc_code_to_slug('HLA-A*02:01:59'), 'hla_a_02_01_59')

    def test_mouse(self):
        self.assertEqual(mhc_code_to_slug('H2-Kb'), 'h2_kb')


class TestMHCSlugToCode(TestCase):
    def test_human(self):
        self.assertEqual(mhc_slug_to_code('hla_a_02_01'), 'HLA-A*02:01')

    def test_mouse(self):
        self.assertEqual(mhc_slug_to_code('h2_kb'), 'H2-Kb')


class TestCentrePad(TestCase):
    def test_same(self):
        self.assertEqual(centre_pad(['x', 'x', 'x', 'x'], 4), ['x', 'x', 'x', 'x'])
        self.assertEqual(centre_pad(['x', 'x', 'x', 'x', 'x'], 5), ['x', 'x', 'x', 'x', 'x'])

    def test_pad_even_matched(self):
        self.assertEqual(centre_pad(['x', 'x', 'x', 'x'], 8), ['-', '-', 'x', 'x', 'x', 'x', '-', '-'])

    def test_pad_even_unmatched(self):
        self.assertEqual(centre_pad(['x', 'x', 'x', 'x'], 7), ['-', '-', 'x', 'x', 'x', 'x', '-'])

    def test_pad_odd_matched(self):
        self.assertEqual(centre_pad(['x', 'x', 'x'], 5), ['-', 'x', 'x', 'x', '-'])

    def test_pad_odd_unmatched(self):
        self.assertEqual(centre_pad(['x', 'x', 'x'], 6), ['-', '-', 'x', 'x', 'x', '-'])

    def test_crop_even_matched(self):
        self.assertEqual(centre_pad(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], 4), ['c', 'd', 'e', 'f'])

    def test_crop_even_unmatched(self):
        self.assertEqual(centre_pad(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], 5), ['c', 'd', 'e', 'f', 'g'])

    def test_crop_odd_matched(self):
        self.assertEqual(centre_pad(['a', 'b', 'c', 'd', 'e', 'f', 'g'], 5), ['b', 'c', 'd', 'e', 'f'])

    def test_crop_odd_unmatched(self):
        self.assertEqual(centre_pad(['a', 'b', 'c', 'd', 'e', 'f', 'g'], 4), ['c', 'd', 'e', 'f'])

    def test_peptide(self):
        self.assertEqual(
            centre_pad(['A', 'P', 'R', 'G', 'P', 'H', 'G', 'G', 'A', 'A', 'S', 'G', 'L'], 12),
            ['P', 'R', 'G', 'P', 'H', 'G', 'G', 'A', 'A', 'S', 'G', 'L'],
        )
