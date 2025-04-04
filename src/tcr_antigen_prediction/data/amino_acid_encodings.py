"""Different encodings for amino acid sequences."""

from collections import OrderedDict

import numpy as np
from Bio.SeqUtils import IUPACData

PROTEIN_LETTERS: list[str] = sorted(IUPACData.protein_letters_3to1.values())
'''Amino acid one letter codes.'''


def _create_one_hot_encoding() -> dict[str, np.ndarray]:
    one_hot_mapping = {}

    zero_arr = np.zeros(len(PROTEIN_LETTERS), dtype=int)

    for idx, olc in enumerate(PROTEIN_LETTERS):
        encoding = zero_arr.copy()
        encoding[idx] = 1
        one_hot_mapping[olc] = encoding

    one_hot_mapping['-'] = zero_arr.copy()

    return one_hot_mapping


ONE_HOT_ENCODING: dict[str, np.ndarray] = _create_one_hot_encoding()
'''One hot encoding of amino acid letters.'''

BLOSUM_50_ENCODING: OrderedDict[str, np.ndarray] = OrderedDict()
'''BLOSUM 50 Encoding of amino acid substitutions based off of DOI: 10.1109/TENCONSpring.2014.6862994.'''

BLOSUM_50_ENCODING['A'] = np.array([5, -1, -2, -1, -3, 0, -2, -1, -1, -2, -1, -1, -1, -1, -2, 1, 0, 0, -3, -2])
BLOSUM_50_ENCODING['C'] = np.array([-1, 13, -4, -3, -2, -3, -3, -2, -3, -2, -2, -2, -4, -3, -4, -1, -1, -1, -5, -3])
BLOSUM_50_ENCODING['D'] = np.array([-2, -4, 8, 2, -5, -1, -1, -4, -1, -4, -4, 2, -1, 0, -2, 0, -1, -4, -5, -3])
BLOSUM_50_ENCODING['E'] = np.array([-1, -3, 2, 6, -3, -3, 0, -4, 1, -3, -2, 0, -1, 2, 0, -1, -1, -3, -3, -2])
BLOSUM_50_ENCODING['F'] = np.array([-3, -2, -5, -3, 8, -4, -1, 0, -4, 1, 0, -4, -4, -4, -3, -3, -2, -1, 1, 4])
BLOSUM_50_ENCODING['G'] = np.array([0, -3, -1, -3, -4, 8, -2, -4, -2, -4, -3, 0, -2, -2, -3, 0, -2, -4, -3, -3])
BLOSUM_50_ENCODING['H'] = np.array([-2, -3, -1, 0, -1, -2, 10, -4, 0, -3, -1, 1, -2, 1, 0, -1, -2, -4, -3, 2])
BLOSUM_50_ENCODING['I'] = np.array([-1, -2, -4, -4, 0, -4, -4, 5, -3, 2, 2, -3, -3, -3, -4, -3, -1, 4, -3, -1])
BLOSUM_50_ENCODING['K'] = np.array([-1, -3, -1, 1, -4, -2, 0, -3, 6, -3, -2, 0, -1, 2, 3, 0, -1, -3, -3, -2])
BLOSUM_50_ENCODING['L'] = np.array([-2, -2, -4, -3, 1, -4, -3, 2, -3, 5, 3, -4, -4, -2, -3, -3, -1, 1, -2, -1])
BLOSUM_50_ENCODING['M'] = np.array([-1, -2, -4, -2, 0, -3, -1, 2, -2, 3, 7, -2, -3, 0, -2, -2, -1, 1, -1, 0])
BLOSUM_50_ENCODING['N'] = np.array([-1, -2, 2, 0, -4, 0, 1, -3, 0, -4, -2, 7, -2, 0, -1, 1, 0, -3, -4, -2])
BLOSUM_50_ENCODING['P'] = np.array([-1, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10, -1, -3, -1, -1, -3, -4, -3])
BLOSUM_50_ENCODING['Q'] = np.array([-1, -3, 0, 2, -4, -2, 1, -3, 2, -2, 0, 0, -1, 7, 1, 0, -1, -3, -1, -1])
BLOSUM_50_ENCODING['R'] = np.array([-2, -4, -2, 0, -3, -3, 0, -4, 3, -3, -2, -1, -3, 1, 7, -1, -1, -3, -3, -1])
BLOSUM_50_ENCODING['S'] = np.array([1, -1, 0, -1, -3, 0, -1, -3, 0, -3, -2, 1, -1, 0, -1, 5, 2, -2, -4, -2])
BLOSUM_50_ENCODING['T'] = np.array([0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 2, 5, 0, -3, -2])
BLOSUM_50_ENCODING['V'] = np.array([0, -1, -4, -3, -1, -4, -4, 4, -3, 1, 1, -3, -3, -3, -3, -2, 0, 5, -3, -1])
BLOSUM_50_ENCODING['W'] = np.array([-3, -5, -5, -3, 1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4, -3, -3, 15, 2])
BLOSUM_50_ENCODING['Y'] = np.array([-2, -3, -3, -2, 4, -3, 2, -1, -2, -1, 0, -2, -3, -1, -1, -2, -2, -1, 2, 8])
BLOSUM_50_ENCODING['-'] = np.array([-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5])
