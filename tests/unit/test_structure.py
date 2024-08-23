import os
from unittest import TestCase

from Bio.PDB import PDBParser

from tcr_antigen_prediction.structure import get_sequence, extract_chains

TEST_DATA = 'tests/data'


class TestGetSequence(TestCase):
    def setUp(self):
        parser = PDBParser()
        self.structure = parser.get_structure('test', os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

    def test(self):
        sequence = get_sequence(self.structure, 'D')
        self.assertEqual(sequence, ('KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSD'
                                    'SATYLCAVTTDSWGKLQFGAGTQVVVTP'))

    def test_res_selection(self):
        sequence = get_sequence(self.structure, 'D', residue_range=set(range(105, 117 + 1)))
        self.assertEqual(sequence, ('AVTTDSWGKLQ'))


class TestExtractChains(TestCase):
    def setUp(self):
        parser = PDBParser()
        self.structure = parser.get_structure('test', os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

    def test_E(self):
        test_struct = extract_chains(self.structure, ['E'])
        test_chains = set([chain.id for chain in test_struct[0]])

        self.assertEqual(test_chains, {'E'})

    def test_DE(self):
        test_struct = extract_chains(self.structure, ['D', 'E'])
        test_chains = set([chain.id for chain in test_struct[0]])

        self.assertEqual(test_chains, {'D', 'E'})

    def test_DECA(self):
        test_struct = extract_chains(self.structure, ['D', 'E', 'C', 'A'])
        test_chains = set([chain.id for chain in test_struct[0]])

        self.assertEqual(test_chains, {'D', 'E', 'C', 'A'})