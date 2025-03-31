from unittest import TestCase

import numpy as np
import pandas as pd

from tcr_antigen_prediction.data.utils import (
    assign_mhc_class,
    assign_species,
    centre_pad,
    create_shared_groups,
    find_common_groups,
    get_cdr_sequences,
    get_mhc_pseudo_sequence,
    left_pad,
    mhc_code_to_slug,
    mhc_slug_to_code,
    right_pad,
    stitch_sequence,
)


class TestStitchSequence(TestCase):
    def test_alpha_chain(self):
        self.assertEqual(
            stitch_sequence('TRAV6N-7*01', 'TRAJ11*01', 'CALAPPDKLTF', 'Mouse'),
            'MNSSPGFMTVMLLIFTRAHGDSVTQTEGQVALSEEDFLTIHCNYSASGYPALFWYVQYPGEGPQFLFRASRDKEKGSSRGFEATYDKGTTSFHLRKASVQESDSAV'
            'YYCALAPPDKLTFGKGTVLLVSPDIQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKTVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNAT'
            'YPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS',
        )

    def test_beta_chain(self):
        self.assertEqual(
            stitch_sequence('TRBV12-1*01', 'TRBJ2-2*01', 'CASSVPGQGDTGQLYF', 'Mouse'),
            'MSNTVLADSAWGITLLSWVTVFLLGTSSADSGVVQSPRHIIKEKGGRSVLTCIPISGHSNVVWYQQTLGKELKFLIQHYEKVERDKGFLPSRFSVQQFDDYHSEMN'
            'MSALELEDSAMYFCASSVPGQGDTGQLYFGEGSKLTVLEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNY'
            'SYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS',
        )


class TestGetCDRSequences(TestCase):
    def test_alpha_chain(self):
        cdr1_alpha, cdr2_alpha, cdr3_alpha = get_cdr_sequences(
            'MDAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASTWGRASTDTQY'
            'FGPGTRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRC'
            'QVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD',
        )

        self.assertEqual(cdr1_alpha, 'SGHNS')
        self.assertEqual(cdr2_alpha, 'FNNNVP')
        self.assertEqual(cdr3_alpha, 'ASTWGRASTDTQY')

    def test_beta_chain(self):
        cdr1_beta, cdr2_beta, cdr3_beta = get_cdr_sequences(
            'MAQTVTQSQPEMSVQEAETVTLSCTYDTSESDYYLFWYKQPPSRQMILVIRQEAYKQQNATENRFSVNFQKAAKSFSLKISDSQLGDAAMYFCASSGNTPLVFGKG'
            'TRLSVIPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS'
        )

        self.assertEqual(cdr1_beta, 'TSESDYY')
        self.assertEqual(cdr2_beta, 'QEAYKQQN')
        self.assertEqual(cdr3_beta, 'ASSGNTPLV')

    def test_invalid(self):
        cdr1, cdr2, cdr3 = get_cdr_sequences(
            'CLFTDFDSQTNVSQSKDSDVNSIIPEDTRLSVIPNSMDFKSNSAVAWSNKSDFACANAIQNPDPAVYQLRDSKSSYITDKCVLDMRFNDKSVTFFPSPESSMAQTV'
            'TQSQPEMLVIRQEAYKQQNATENRFSMYFCSVQEAETVTLSCTYDTSVNFQKAAKSFSLKISDSQLGDAAESDYYLFWYKQPPSRQMIASSGNTPLVFGKG'
        )

        self.assertIsNone(cdr1)
        self.assertIsNone(cdr2)
        self.assertIsNone(cdr3)


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


class TestAssignMHCClass(TestCase):
    def test_class_i_human(self):
        self.assertEqual(assign_mhc_class('HLA-A*02:01'), 'MH1')

    def test_class_ii_human(self):
        self.assertEqual(assign_mhc_class('HLA-DQ*02:01'), 'MH2')

    def test_class_i_mouse(self):
        self.assertEqual(assign_mhc_class('H2-D'), 'MH1')

    def test_class_ii_mouse(self):
        self.assertEqual(assign_mhc_class('H2-IA'), 'MH2')


class TestGetMHCPseudoSequence(TestCase):
    def setUp(self):
        self.imgt_pseudo_seq_positions = {
            'alpha': ['62', '63', '65', '66', '68', '69', '70', '72', '73', '75', '76', '79'],
            'beta': ['58', '61A', '62', '63', '65', '66', '67', '69', '70', '72', '72A', '73', '76', '77'],
        }

    def test_mh1(self):
        mhc_sequence = (
            'GSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPD'
            'GRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAE'
            'ITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE'
        )

        self.assertEqual(
            get_mhc_pseudo_sequence(mhc_sequence, None, 'MH1', self.imgt_pseudo_seq_positions),
            'RNQIKAQQTRERKAREEQRAYEGEEW',
        )

    def test_mh2(self):
        mhc1_sequence = (
            'IKEEHTIIQAEFYLLPDKRGEFMFDFDGDEIFHVDIEKSETIWRLEEFAKFASFEAQGALANIAVDKANLDVMKERSNNTPDANVAPEVTVLSRSPVNLGEPNILI'
            'CFIDKFSPPVVNVTWLRNGRPVTEGVSETVFLPRDDHLFRKFHYLTFLPSTDDFYDCEVDHWGLEEPLRKTWEFEEKTLLPETKEN'
        )

        mhc2_sequence = (
            'RDSRDRMVNHFIAEFKRKGGSLVPRGSGGGGSRPWFLEYCKSECHFYNGTQRVRLLVRYFYNLEENLRFDSDVGEFRAVTELGRPDAENWNSQPEFLEQKRAEVDT'
            'VCRHNYEIFDNFLVPRRVEPTVTVYPTKTQPLEHHNLLVCSVSDFYPGNIEVRWFRNGKEEKTGIVSTGLVRNGDWTFQTLVMLETVPQSGEVYTCQVEHPSLTDP'
            'VTVEWKAQSTSAQNK'
        )

        self.assertEqual(
            get_mhc_pseudo_sequence(mhc1_sequence, mhc2_sequence, 'MH2', self.imgt_pseudo_seq_positions),
            'FEQGLANAVKADNQEFEQKAEDTVHN',
        )


class TestAssignSpecies(TestCase):
    def test_human(self):
        self.assertEqual(assign_species('HLA-A*02:01'), 'Human')
        self.assertEqual(assign_species('HLA-E*03:01'), 'Human')
        self.assertEqual(assign_species('HLA-DQ*02:01'), 'Human')

    def test_mouse(self):
        self.assertEqual(assign_species('H2-A'), 'Mouse')
        self.assertEqual(assign_species('H2-D'), 'Mouse')

    def test_chicken(self):
        self.assertEqual(assign_species('Gaga-BF1*002:01:01'), 'Chicken')


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


class TestRightPad(TestCase):
    def test_pad(self):
        self.assertEqual(right_pad(['A', 'B', 'C'], 5), ['A', 'B', 'C', '-', '-'])

    def test_crop(self):
        self.assertEqual(right_pad(['A', 'B', 'C'], 2), ['A', 'B'])

    def test_nothing(self):
        self.assertEqual(right_pad(['A', 'B', 'C'], 3), ['A', 'B', 'C'])


class TestLeftPad(TestCase):
    def test_pad(self):
        self.assertEqual(left_pad(['A', 'B', 'C'], 5), ['-', '-', 'A', 'B', 'C'])

    def test_crop(self):
        self.assertEqual(left_pad(['A', 'B', 'C'], 2), ['B', 'C'])

    def test_nothing(self):
        self.assertEqual(left_pad(['A', 'B', 'C'], 3), ['A', 'B', 'C'])


class TestFindCommonGroups(TestCase):
    def test_none(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        with self.assertRaises(ValueError):
            find_common_groups(data, [])

    def test_1d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 2, 3, 3, 3, 4, 4, 5, 5]))

    def test_2d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_b': ['l', 'm', 'l', 'o', 'p', 'q', 'q', 'r', 's', 't'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a', 'col_b'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 1, 2, 2, 2, 2, 2, 3, 3]))

    def test_4d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_b': ['l', 'm', 'l', 'o', 'p', 'q', 'q', 'r', 's', 't'],
                'col_c': ['x', 'y', 'z', 'aa', 'bb', 'cc', 'dd', 'ee', 'ee', 'ee'],
                'col_d': ['aaa', 'bbb', 'ccc', 'ddd', 'eee', 'fff', 'ggg', 'ggg', 'ggg', 'ggg'],
                'col_e': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a', 'col_b', 'col_c', 'col_d'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 1, 2, 2, 2, 2, 2, 2, 2]))


class TestCreateSharedGroups(TestCase):
    def test(self):
        sequences = pd.DataFrame(
            {
                'peptide_sequence': ['a', 'b', 'a', 'c', 'd'],
                'collated_cdrs': ['A', 'A', 'B', 'C', 'C'],
                'dummy': [1, 2, 4, 3, 1],
            }
        )

        groups = create_shared_groups(sequences)

        self.assertEqual(len(groups), 2)
        pd.testing.assert_frame_equal(
            groups[0],
            pd.DataFrame(
                {
                    'peptide_sequence': ['a', 'b', 'a'],
                    'collated_cdrs': ['A', 'A', 'B'],
                    'dummy': [1, 2, 4],
                }
            ),
        )
        pd.testing.assert_frame_equal(
            groups[1],
            pd.DataFrame(
                {
                    'peptide_sequence': ['c', 'd'],
                    'collated_cdrs': ['C', 'C'],
                    'dummy': [3, 1],
                },
                index=[3, 4],
            ),
        )
