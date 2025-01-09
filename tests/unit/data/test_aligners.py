from unittest import TestCase

from tcr_antigen_prediction.data.aligners import align_sequences


class TestAlignSequence(TestCase):
    def test(self):
        alignment, score = align_sequences('GCATGCG', 'GATTACA')

        self.assertEqual(
            alignment, [('G', 'G'), ('C', '-'), ('A', 'A'), ('T', 'T'), ('G', 'T'), ('-', 'A'), ('C', 'C'), ('G', 'A')]
        )

        self.assertEqual(score, 0.0)
