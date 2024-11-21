import os
from unittest import TestCase

from Bio.PDB import PDBParser

from tcr_antigen_prediction.imgt_numbering import IMGT_MH2_ABD, IMGT_VARIABLE_DOMAIN
from tcr_antigen_prediction.missing_residues import (
    get_missing_atoms,
    get_missing_residues,
    screen_chain,
    trim_start_and_end,
)
from tcr_antigen_prediction.structure import get_header

TEST_DATA = 'tests/data'


class TestGetMissingAtoms(TestCase):
    def test_3qiw(self):
        with open(os.path.join(TEST_DATA, '3qiw_raw.pdb')) as fh:
            header = get_header(fh.read())

        missing_atoms = get_missing_atoms(header)

        self.assertEqual(
            missing_atoms, [{'residue_name': 'LYS', 'chain_id': 'A', 'residue_seq_id': 38, 'residue_insert_code': ''}]
        )


class TestGetMissingResidues(TestCase):
    def test_3qiw(self):
        with open(os.path.join(TEST_DATA, '3qiw_raw.pdb')) as fh:
            header = get_header(fh.read())

        missing_residues = get_missing_residues(header)
        self.assertEqual(
            missing_residues,
            [
                {'residue_name': 'ILE', 'chain_id': 'A', 'residue_seq_id': 1, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'A', 'residue_seq_id': 2, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'A', 'residue_seq_id': 182, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'A', 'residue_seq_id': 183, 'residue_insert_code': ''},
                {'residue_name': 'THR', 'chain_id': 'A', 'residue_seq_id': 184, 'residue_insert_code': ''},
                {'residue_name': 'LEU', 'chain_id': 'A', 'residue_seq_id': 185, 'residue_insert_code': ''},
                {'residue_name': 'LEU', 'chain_id': 'A', 'residue_seq_id': 186, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'A', 'residue_seq_id': 187, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'A', 'residue_seq_id': 188, 'residue_insert_code': ''},
                {'residue_name': 'THR', 'chain_id': 'A', 'residue_seq_id': 189, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'A', 'residue_seq_id': 190, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'A', 'residue_seq_id': 191, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'A', 'residue_seq_id': 192, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'B', 'residue_seq_id': 105, 'residue_insert_code': ''},
                {'residue_name': 'THR', 'chain_id': 'B', 'residue_seq_id': 106, 'residue_insert_code': ''},
                {'residue_name': 'GLN', 'chain_id': 'B', 'residue_seq_id': 107, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'B', 'residue_seq_id': 108, 'residue_insert_code': ''},
                {'residue_name': 'LEU', 'chain_id': 'B', 'residue_seq_id': 109, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'B', 'residue_seq_id': 110, 'residue_insert_code': ''},
                {'residue_name': 'HIS', 'chain_id': 'B', 'residue_seq_id': 111, 'residue_insert_code': ''},
                {'residue_name': 'HIS', 'chain_id': 'B', 'residue_seq_id': 112, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'B', 'residue_seq_id': 113, 'residue_insert_code': ''},
                {'residue_name': 'ARG', 'chain_id': 'B', 'residue_seq_id': 133, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'B', 'residue_seq_id': 134, 'residue_insert_code': ''},
                {'residue_name': 'GLY', 'chain_id': 'B', 'residue_seq_id': 135, 'residue_insert_code': ''},
                {'residue_name': 'VAL', 'chain_id': 'B', 'residue_seq_id': 164, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'B', 'residue_seq_id': 165, 'residue_insert_code': ''},
                {'residue_name': 'GLN', 'chain_id': 'B', 'residue_seq_id': 166, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'B', 'residue_seq_id': 167, 'residue_insert_code': ''},
                {'residue_name': 'GLY', 'chain_id': 'B', 'residue_seq_id': 168, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'B', 'residue_seq_id': 169, 'residue_insert_code': ''},
                {'residue_name': 'VAL', 'chain_id': 'B', 'residue_seq_id': 170, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'B', 'residue_seq_id': 187, 'residue_insert_code': ''},
                {'residue_name': 'TRP', 'chain_id': 'B', 'residue_seq_id': 188, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'B', 'residue_seq_id': 189, 'residue_insert_code': ''},
                {'residue_name': 'ALA', 'chain_id': 'B', 'residue_seq_id': 190, 'residue_insert_code': ''},
                {'residue_name': 'GLN', 'chain_id': 'B', 'residue_seq_id': 191, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'B', 'residue_seq_id': 192, 'residue_insert_code': ''},
                {'residue_name': 'THR', 'chain_id': 'B', 'residue_seq_id': 193, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'B', 'residue_seq_id': 194, 'residue_insert_code': ''},
                {'residue_name': 'ALA', 'chain_id': 'B', 'residue_seq_id': 195, 'residue_insert_code': ''},
                {'residue_name': 'GLN', 'chain_id': 'B', 'residue_seq_id': 196, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'B', 'residue_seq_id': 197, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'B', 'residue_seq_id': 198, 'residue_insert_code': ''},
                {'residue_name': 'GLY', 'chain_id': 'C', 'residue_seq_id': 0, 'residue_insert_code': ''},
                {'residue_name': 'ASP', 'chain_id': 'C', 'residue_seq_id': 1, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 145, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'C', 'residue_seq_id': 146, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'C', 'residue_seq_id': 175, 'residue_insert_code': ''},
                {'residue_name': 'LYS', 'chain_id': 'C', 'residue_seq_id': 176, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 177, 'residue_insert_code': ''},
                {'residue_name': 'ASP', 'chain_id': 'C', 'residue_seq_id': 178, 'residue_insert_code': ''},
                {'residue_name': 'PHE', 'chain_id': 'C', 'residue_seq_id': 179, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'C', 'residue_seq_id': 187, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 188, 'residue_insert_code': ''},
                {'residue_name': 'ILE', 'chain_id': 'C', 'residue_seq_id': 189, 'residue_insert_code': ''},
                {'residue_name': 'ILE', 'chain_id': 'C', 'residue_seq_id': 190, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'C', 'residue_seq_id': 191, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'C', 'residue_seq_id': 192, 'residue_insert_code': ''},
                {'residue_name': 'ASP', 'chain_id': 'C', 'residue_seq_id': 193, 'residue_insert_code': ''},
                {'residue_name': 'THR', 'chain_id': 'C', 'residue_seq_id': 194, 'residue_insert_code': ''},
                {'residue_name': 'PHE', 'chain_id': 'C', 'residue_seq_id': 195, 'residue_insert_code': ''},
                {'residue_name': 'PHE', 'chain_id': 'C', 'residue_seq_id': 196, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'C', 'residue_seq_id': 197, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 198, 'residue_insert_code': ''},
                {'residue_name': 'PRO', 'chain_id': 'C', 'residue_seq_id': 199, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'C', 'residue_seq_id': 200, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 201, 'residue_insert_code': ''},
                {'residue_name': 'SER', 'chain_id': 'C', 'residue_seq_id': 202, 'residue_insert_code': ''},
                {'residue_name': 'LEU', 'chain_id': 'D', 'residue_seq_id': 183, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'D', 'residue_seq_id': 184, 'residue_insert_code': ''},
                {'residue_name': 'GLU', 'chain_id': 'D', 'residue_seq_id': 219, 'residue_insert_code': ''},
                {'residue_name': 'ASN', 'chain_id': 'D', 'residue_seq_id': 220, 'residue_insert_code': ''},
                {'residue_name': 'ASP', 'chain_id': 'D', 'residue_seq_id': 244, 'residue_insert_code': ''},
            ],
        )


class TestTrimStartAndEnd(TestCase):
    def test(self):
        sequence = [
            ['A', 1, '', 'missing'],
            ['A', 2, '', 'missing'],
            ['A', 3, '', 'not-missing'],
            ['A', 4, '', 'not-missing'],
            ['A', 5, '', 'not-missing'],
            ['A', 5, 'A', 'missing'],
            ['A', 6, '', 'not-missing'],
            ['A', 7, '', 'not-missing'],
            ['A', 8, '', 'not-missing'],
            ['A', 9, '', 'missing'],
            ['A', 10, '', 'missing'],
        ]

        output_seq = trim_start_and_end(sequence)

        self.assertEqual(
            output_seq,
            [
                ['A', 3, '', 'not-missing'],
                ['A', 4, '', 'not-missing'],
                ['A', 5, '', 'not-missing'],
                ['A', 5, 'A', 'missing'],
                ['A', 6, '', 'not-missing'],
                ['A', 7, '', 'not-missing'],
                ['A', 8, '', 'not-missing'],
            ],
        )


class TestScreenChain(TestCase):
    def test_3qiw_pass(self):
        with open(os.path.join(TEST_DATA, '3qiw_raw.pdb')) as fh:
            header = get_header(fh.read())

        missing_residues = get_missing_residues(header)
        missing_atoms = get_missing_atoms(header)

        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('imgt', os.path.join(TEST_DATA, '3qiw_imgt.pdb'))
        raw_structure = pdb_parser.get_structure('raw', os.path.join(TEST_DATA, '3qiw_raw.pdb'))

        self.assertTrue(
            screen_chain(structure, raw_structure, 'C', missing_residues, missing_atoms, IMGT_VARIABLE_DOMAIN)
        )

    def test_3qiw_fail(self):
        with open(os.path.join(TEST_DATA, '3qiw_raw.pdb')) as fh:
            header = get_header(fh.read())

        missing_residues = get_missing_residues(header)
        missing_atoms = get_missing_atoms(header)

        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('imgt', os.path.join(TEST_DATA, '3qiw_imgt.pdb'))
        raw_structure = pdb_parser.get_structure('raw', os.path.join(TEST_DATA, '3qiw_raw.pdb'))

        self.assertFalse(screen_chain(structure, raw_structure, 'A', missing_residues, missing_atoms, IMGT_MH2_ABD))
