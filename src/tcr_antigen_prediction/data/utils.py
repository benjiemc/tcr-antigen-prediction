"""Utility functions for processing TCR:pMHC data."""

import logging
import random
import re
import typing
import warnings

import pandas as pd
from Stitchr import stitchr as st
from Stitchr import stitchrfunctions as fxn

from tcr_antigen_prediction.data.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3, MHC_I_IMGT_BETA_HELIX_START

logger = logging.getLogger(__name__)

try:
    import anarci
except ImportError:
    logger.exception(
        'Some functions require ANARCI which is unavailable. see here for install instructions '
        'https://github.com/oxpig/ANARCI.'
    )


def mhc_code_to_slug(code: str) -> str:
    """Convert MHC allele codes into slugs.

    >>> mhc_code_to_slug('HLA-A*02:01:59')
    'hla_a_02_01_59'

    """
    slug = code.lower()
    slug = re.sub(r'[*:-]', '_', slug)

    return slug


def mhc_slug_to_code(slug: str) -> str:
    """Convert mhc slugs into the allele codes.

    >>> mhc_slug_to_code('hla_a_02_01')
    'HLA-A*02:01'

    >>> mhc_slug_to_code('h2_kb')
    'H2-Kb'

    TODO: Make this better for other mouse alleles
    """
    species = 'human' if slug.startswith('hla') else 'mouse' if slug.startswith('h2') else None

    if species == 'human':
        code = slug.upper()
        code = code.split('_')
        code = code[0] + '-' + code[1] + '*' + ':'.join(code[2:])

    elif species == 'mouse':
        code = slug.title()
        code = code.replace('_', '-')

    else:
        logger.error('Species not found, outputting slug')
        code = slug

    return code


def assign_mhc_class(allele: str) -> str | None:
    """Assign MHC Class based on the allele code provided."""
    if re.search('HLA-[A-CE]', allele) or re.search('H2-[DKLbdkq]', allele):
        return 'MH1'

    if re.search('HLA-D[PMOQR][AB]?', allele) or re.search('H2-I[AE]', allele):
        return 'MH2'

    return None


def assign_species(allele: str) -> str | None:
    """Assign species based on MHC allele code."""
    if re.search('^HLA-', allele):
        return 'Human'

    if re.search('^H2-?', allele):
        return 'Mouse'

    if re.search('^Gaga', allele):
        return 'Chicken'

    return None


def stitch_sequence(v_gene: str, j_gene: str, cdr3: str, species: str) -> str | None:
    """Create full length TCR sequence from V gene, J, gene, CDR3 and species information."""

    def create_input_args(args: dict, gene_types: list) -> tuple:
        input_args, chain = fxn.sort_input(args)
        codons = fxn.get_optimal_codons(input_args['codon_usage_path'], input_args['species'])
        imgt_dat, tcr_functionality, partial = fxn.get_imgt_data(chain, gene_types, input_args['species'])

        if input_args['extra_genes']:
            imgt_dat, tcr_functionality = fxn.get_additional_genes(imgt_dat, tcr_functionality)
            input_args['skip_c_checks'] = True

        if input_args['preferred_alleles_path']:
            preferred_alleles = fxn.get_preferred_alleles(
                input_args['preferred_alleles_path'],
                gene_types,
                imgt_dat,
                partial,
                chain,
            )

        else:
            preferred_alleles = {}

        return input_args, imgt_dat, tcr_functionality, partial, codons, preferred_alleles

    gene_types = list(fxn.regions.values())

    args = {
        'v': v_gene,
        'j': j_gene,
        'cdr3': cdr3,
        'species': species.upper(),
        'c': '',
        'l': '',
        'aa': '',
        'name': '',
        'seamless': False,
        '5_prime_seq': '',
        '3_prime_seq': '',
        'extra_genes': False,
        'preferred_alleles_path': '',
        'codon_usage_path': '',
        'j_warning_threshold': 3,
        'skip_c_checks': False,
        'suppress_warnings': False,
    }

    try:
        (
            input_args,
            imgt_dat,
            tcr_functionality,
            partial,
            codons,
            preferred_alleles,
        ) = create_input_args(args, gene_types)

    except ValueError as err:
        logger.debug(err, args)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _, stitched, offset = st.stitch(
                input_args,
                imgt_dat,
                tcr_functionality,
                partial,
                codons,
                input_args['j_warning_threshold'],
                preferred_alleles,
            )

    except ValueError:
        if cdr3[-1] == 'F':
            cdr3 = cdr3[:-1] + 'W'

        elif cdr3[-1] == 'W':
            cdr3 = cdr3[:-1] + 'F'

        args['cdr3'] = cdr3

        (
            input_args,
            imgt_dat,
            tcr_functionality,
            partial,
            codons,
            preferred_alleles,
        ) = create_input_args(args, gene_types)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _, stitched, offset = st.stitch(
                    input_args,
                    imgt_dat,
                    tcr_functionality,
                    partial,
                    codons,
                    input_args['j_warning_threshold'],
                    preferred_alleles,
                )

        except ValueError:
            return None

    except Exception:  # noqa: BLE001 Stitchr will throw a base exception sometimes
        return None

    return fxn.translate_nt('N' * offset + stitched)


def get_cdr_sequences(sequence: str) -> tuple[str | None, str | None, str | None]:
    """Get CDR sequences using Anarci."""
    try:
        numbering, _ = anarci.number(sequence)

    except AssertionError as err:
        logger.debug(err, sequence)
        return None, None, None

    if not numbering:
        logger.debug('Sequence not numberable: %s', sequence)
        return None, None, None

    if len(numbering) == 2:  # noqa: PLR2004 Sometimes anarci will find two valid numberings
        numbering = numbering[0]

    cdr1 = ''.join([letter for (seq_id, _), letter in numbering if letter != '-' and seq_id in IMGT_CDR1])
    cdr2 = ''.join([letter for (seq_id, _), letter in numbering if letter != '-' and seq_id in IMGT_CDR2])
    cdr3 = ''.join([letter for (seq_id, _), letter in numbering if letter != '-' and seq_id in IMGT_CDR3])

    return cdr1, cdr2, cdr3


def get_mhc_pseudo_sequence(  # noqa: C901, PLR0912 let this function be complicated
    mhc1_sequence: str,
    mhc2_sequence: str | None,
    mhc_type: str,
    mhc_pseudo_seq_imgt_positions: dict[str, list[str]],
) -> list[str]:
    """Get the MHC pseudo sequence using anarci and imgt mhc pseudo sequence positiions."""
    mhc_pseudo_seq_imgt_positions = {
        helix: [
            (
                int(''.join([char for char in resi if char.isnumeric()])),
                ''.join([char for char in resi if not char.isnumeric()]),
            )
            for resi in residues
        ]
        for helix, residues in mhc_pseudo_seq_imgt_positions.items()
    }
    numberings = {}

    numbering1, _ = anarci.number(mhc1_sequence)
    numberings['chain1'] = numbering1

    if not pd.isna(mhc2_sequence):
        numbering2, _ = anarci.number(mhc2_sequence)
        numberings['chain2'] = numbering2

    mhc_pseudo_seq = []
    for helix, residues in mhc_pseudo_seq_imgt_positions.items():
        for imgt_seq_id, imgt_insert_code in residues:
            if helix == 'alpha':
                for (seq_id, insert_code), res_olc in numberings['chain1']:
                    if seq_id == imgt_seq_id and insert_code.strip() == imgt_insert_code:
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

            elif helix == 'beta' and mhc_type == 'MH1':
                for (seq_id, insert_code), res_olc in numberings['chain1']:
                    if (
                        seq_id == (imgt_seq_id + MHC_I_IMGT_BETA_HELIX_START)
                        and insert_code.strip() == imgt_insert_code
                    ):
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

            elif helix == 'beta' and mhc_type == 'MH2':
                for (seq_id, insert_code), res_olc in numberings['chain2']:
                    if seq_id == imgt_seq_id and insert_code.strip() == imgt_insert_code:
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

    return ''.join(mhc_pseudo_seq)


def centre_pad(sequence: list[str], pad_length: int) -> list[str]:
    """Pad or crop sequence to pad_length.

    Centre pad a sequence with '-'s to a specified pad length or crop sequences to size if they are longer than the
    specified pad_length.

    """
    sequence_length = len(sequence)

    if sequence_length < pad_length:
        if (pad_length % 2) == (sequence_length % 2):
            pad_per_side = (pad_length - sequence_length) // 2
            pad_right_side = pad_per_side
            pad_left_side = pad_per_side

        else:
            pad_right_side = (pad_length - sequence_length) // 2
            pad_left_side = pad_right_side + 1

        return (['-'] * pad_left_side) + sequence + (['-'] * pad_right_side)

    if sequence_length > pad_length:
        if (pad_length % 2) == (sequence_length % 2):
            crop_per_side = (sequence_length - pad_length) // 2
            crop_right_side = crop_per_side
            crop_left_side = crop_per_side

        else:
            crop_right_side = (sequence_length - pad_length) // 2
            crop_left_side = crop_right_side + 1

        return sequence[crop_left_side:-crop_right_side] if crop_right_side > 0 else sequence[crop_left_side:]

    return sequence


def right_pad(sequence: list[str], pad_length: int) -> list[str]:
    """Right pad or crop sequence based on the specified pad length.

    Args:
        sequence: sequence to pad
        pad_length: length to pad or crop to

    Returns:
        sequence with '-' used to represent null pad values

    Examples:
        >>> right_pad(['A', 'B', 'C'], 5)
        ['A', 'B', 'C', '-', '-']

        >>> right_pad(['A', 'B', 'C'], 2)
        ['A', 'B']

        >>> right_pad(['A', 'B', 'C'], 3)
        ['A', 'B', 'C']

    """
    sequence_length = len(sequence)

    if sequence_length < pad_length:
        num_to_add = pad_length - sequence_length
        return sequence + ['-'] * num_to_add

    if sequence_length > pad_length:
        num_to_remove = sequence_length - pad_length
        return sequence[:-num_to_remove]

    return sequence


def left_pad(sequence: list[str], pad_length: int) -> list[str]:
    """Left pad or crop sequence based on the specified pad length.

    Args:
        sequence: sequence to pad
        pad_length: length to pad or crop to

    Returns:
        sequence with '-' used to represent null pad values

    Examples:
        >>> left_pad(['A', 'B', 'C'], 5)
        ['-', '-', 'A', 'B', 'C']

        >>> left_pad(['A', 'B', 'C'], 2)
        ['B', 'C']

        >>> left_pad(['A', 'B', 'C'], 3)
        ['A', 'B', 'C']

    """
    sequence_length = len(sequence)

    if sequence_length < pad_length:
        num_to_add = pad_length - sequence_length
        return ['-'] * num_to_add + sequence

    if sequence_length > pad_length:
        num_to_remove = sequence_length - pad_length
        return sequence[num_to_remove:]

    return sequence


def create_even_folds(
    counts: list[tuple[typing.Any, int]],
    num_folds: int = 5,
    seed: int | None = None,
) -> tuple[list[int]]:
    """Create roughly even data folds based on the number of data points belonging to a group.

    Args:
        counts: list containing tuples with the group object and the count of obejcts belonging to that group
        num_folds: number of folds to make in the dataset
        seed: optional seed for random shuffling of the data

    Returns:
        tuple with the group objects placed in each fold

    """
    if seed:
        random.seed(seed)

    total_sum = sum([size for _, size in counts])
    target_sum = total_sum / num_folds

    groups = [[] for _ in range(num_folds)]
    sum_groups = [0] * num_folds

    random.shuffle(counts)

    for group, size in counts:
        # Find the group with the smallest current sum and add the number
        min_sum_index = min(range(num_folds), key=lambda i: sum_groups[i])

        if sum_groups[min_sum_index] + size <= target_sum:
            groups[min_sum_index].append(group)
            sum_groups[min_sum_index] += size

        else:
            # If adding the number exceeds the target sum, add to the next group
            for i in range(num_folds):
                if i != min_sum_index and sum_groups[i] + size <= target_sum:
                    groups[i].append(group)
                    sum_groups[i] += size
                    break
            else:  # If it doesn't fit anywhere, add it to the original trial
                groups[min_sum_index].append(group)
                sum_groups[min_sum_index] += size

    return tuple(groups)
