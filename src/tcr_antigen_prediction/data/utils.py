"""Utility functions for processing TCR:pMHC data."""

import random
import re
import typing


def mhc_code_to_slug(code: str) -> str:
    """Convert MHC allele codes into slugs.

    >>> mhc_code_to_slug('HLA-A*02:01:59')
    'hla_a_02_01_59'

    """
    slug = code.lower()
    slug = re.sub(r'[*:-]', '_', slug)

    return slug


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
