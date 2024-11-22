"""Select TCR-pMHC structures from STCRDab."""

import argparse
import logging
import os
import random
import sys
import tempfile

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select
from sklearn.cluster import AgglomerativeClustering

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.comparisons import compute_structural_distances
from tcr_antigen_prediction.imgt_numbering import (
    IMGT_CDR,
    IMGT_CDR1,
    IMGT_CDR2,
    IMGT_CDR3,
    IMGT_FRAMEWORK_REGION,
    IMGT_MH1_ABD,
    IMGT_MH2_ABD,
    renumber_chain,
)
from tcr_antigen_prediction.missing_residues import get_alignment, get_missing_atoms, get_missing_residues, screen_chain
from tcr_antigen_prediction.missing_residues.fix import predict_missing_residues
from tcr_antigen_prediction.structure import crop_structure, extract_chains, get_header, get_sequence, replace_chain

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('stcrdab', help='path to the STCRDab')

parser.add_argument('--output', '-o', help='output path')
parser.add_argument('--seed', default=None, type=int, help='random seed for data splitting')
parser.add_argument(
    '--mhc-class-I-tcr-contact-residues',
    nargs='+',
    default=[],
    help='list of IMGT residue codes that are in contact positions on the MHC class I molecules',
)
parser.add_argument(
    '--mhc-class-II-alpha-chain-tcr-contact-residues',
    nargs='+',
    default=[],
    help='list of IMGT residue codes that are in contact positions on the MHC class II alpha-chain',
)
parser.add_argument(
    '--mhc-class-II-beta-chain-tcr-contact-residues',
    nargs='+',
    default=[],
    help='list of IMGT residue codes that are in contact positions on the MHC class II beta-chain',
)

data_split_group = parser.add_argument_group('Data Splits')
data_split_group.add_argument(
    '--train-split', type=float, default=0.70, help='proportion of data to assign to training (Default: 0.70)'
)
data_split_group.add_argument(
    '--validation-split', type=float, default=0.15, help='proportion of data to assign to validation (Default: 0.15)'
)
data_split_group.add_argument(
    '--test-split', type=float, default=0.15, help='proportion of data to assign to testing (Default: 0.15)'
)
data_split_group.add_argument(
    '--pdb-ids-to-exclude', nargs='+', default=None, help='PDB IDs to exclude from the validation/testing data splits'
)

structure_type_group = parser.add_argument_group('Structure Types')
structure_type_group.add_argument(
    '--tcr-types',
    nargs='+',
    default=['abTCR'],
    help='TCR types allowed in dataset (abTCR and/or gdTCR) (Default: abTCR)',
)
structure_type_group.add_argument(
    '--mhc-types',
    nargs='+',
    default=['MH1'],
    help='MHC types allowed in dataset (CD1, GA, GB, MH1, MH2, MR1) (Default: MH1)',
)
structure_type_group.add_argument(
    '--antigen-types',
    nargs='+',
    default=['peptide'],
    help=('MHC types allowed in dataset (carbohydrate, Hapten, peptide, protein, etc) ' '(Default: peptide)'),
)
structure_type_group.add_argument(
    '--crop-structures',
    action='store_true',
    help=('Crop TCR:pMHC structures to the TCR variable domain and MHC antigen binding ' 'domain.'),
)

quality_group = parser.add_argument_group('Quality Selection')
quality_group.add_argument(
    '--resolution-cutoff', type=float, default=3.50, help='maximum resolution allowed (Default: 3.50)'
)
quality_group.add_argument(
    '--remove-structures-missing-residues',
    action='store_true',
    help=(
        'Remove TCR:pMHC structures with missing residues in the TCR variable region or '
        'pMHC antigen binding domain (including peptide)'
    ),
)
quality_group.add_argument(
    '--fix-structures-missing-residues',
    action='store_true',
    help=(
        'Fix TCR:pMHC structures with missing residues in the TCR Fw region or '
        'MHC antigen binding domain rather than discarding (to be used with '
        '`--remove-structures-missing-residues)'
    ),
)
quality_group.add_argument(
    '--structural-similarity-cutoff',
    type=float,
    default=None,
    help='RMSD threshold for structures with the same CDR and peptide sequences (Default: None)',
)

add_logging_arguments(parser)


class SelectChains(Select):
    """Select chains to output."""

    def __init__(self, *chain_ids):
        self.selected_chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.selected_chain_ids


def screen_for_missing_residues(
    df: pd.DataFrame, *, fix_structures: bool = False, fix_dir: str | None = None
) -> pd.DataFrame:
    """Remove entries missing residues in the TCR variable region or pMHC antigen binding domain."""

    def check_structure(  # noqa: PLR0911 TODO refactor this function
        pdb_id,
        alpha_chain_id,
        beta_chain_id,
        antigen_chain_id,
        mhc_chain1_id,
        mhc_chain2_id,
        mhc_type,
        imgt_file_path,
        raw_file_path,
    ):
        if mhc_type == 'MH1':
            mhc_chains = (mhc_chain1_id,)

        elif mhc_type == 'MH2':
            mhc_chains = (mhc_chain1_id, mhc_chain2_id)

        else:
            msg = f"Invalid MHC type {mhc_type}. MHC type must be 'MH1' or 'MH2'"
            raise ValueError(msg)

        logger.debug(
            'Checking %s, chains: %s, and MHC type %s',
            pdb_id,
            '-'.join([alpha_chain_id, beta_chain_id, antigen_chain_id, *mhc_chains]),
            mhc_type,
        )

        with open(raw_file_path, 'r') as fh:
            header = get_header(fh.read())

        pdb_parser = PDBParser(QUIET=True)
        raw_structure = pdb_parser.get_structure(pdb_id, raw_file_path)
        structure = pdb_parser.get_structure(pdb_id + '_imgt_numbered', imgt_file_path)

        missing_atoms = get_missing_atoms(header)
        missing_residues = get_missing_residues(header)

        fixed = False

        try:
            logger.debug('Checking for missing residues or atoms in antigen chain (chain %s)', antigen_chain_id)
            if antigen_chain_id in (
                {res['chain_id'] for res in missing_residues} | {res['chain_id'] for res in missing_atoms}
            ):
                logger.debug('Missing elements found')
                return 'missing'

            for tcr_chain_id in alpha_chain_id, beta_chain_id:
                logger.debug('Checking for missing residues or atoms in CDRs of chain %s', tcr_chain_id)
                if not screen_chain(
                    structure, raw_structure, tcr_chain_id, missing_residues, missing_atoms, IMGT_CDR, trim_ends=False
                ):
                    logger.debug('Missing elements found')
                    return 'missing'

                logger.debug('Checking for missing residues or atoms in framework region of chain %s', tcr_chain_id)
                if not screen_chain(
                    structure, raw_structure, tcr_chain_id, missing_residues, missing_atoms, IMGT_FRAMEWORK_REGION
                ):
                    logger.debug('Missing elements found')

                    if fix_structures:
                        logger.debug('Fixing chain')
                        alignment = get_alignment(structure, raw_structure, tcr_chain_id, missing_residues)

                        new_chain = predict_missing_residues(alignment, structure, tcr_chain_id)
                        new_chain = renumber_chain(new_chain)
                        structure = replace_chain(structure, new_chain)

                        fixed = True

                    else:
                        return 'missing'

            for mhc_chain_id in mhc_chains:
                logger.debug(
                    'Checking for missing residues or atoms in antigen binding domain of chain %s', mhc_chain_id
                )
                if not screen_chain(
                    structure,
                    raw_structure,
                    mhc_chain_id,
                    missing_residues,
                    missing_atoms,
                    IMGT_MH1_ABD if mhc_type == 'MH1' else IMGT_MH2_ABD,
                ):
                    logger.debug('Missing elements found')

                    if fix_structures:
                        logger.debug('Fixing chain')
                        alignment = get_alignment(structure, raw_structure, mhc_chain_id, missing_residues)

                        new_chain = predict_missing_residues(alignment, structure, mhc_chain_id)
                        new_chain = renumber_chain(new_chain)
                        structure = replace_chain(structure, new_chain)

                        io = PDBIO()
                        io.set_structure(structure)
                        io.save(os.path.join(fix_dir, os.path.basename(imgt_file_path)))

                        fixed = True

                    else:
                        return 'missing'

            if fixed:
                io = PDBIO()
                io.set_structure(structure)
                io.save(os.path.join(fix_dir, os.path.basename(imgt_file_path)))

                return 'fixed'

        except KeyError:
            logger.warning(
                'Chain ID not found in raw structure of %s and chains %s. Structure will be removed.',
                pdb_id,
                '-'.join([alpha_chain_id, beta_chain_id, antigen_chain_id, *mhc_chains]),
            )
            return 'unknown'

        else:
            return 'complete'

    structure_status = df.apply(
        lambda row: check_structure(
            row.pdb,
            row.Achain,
            row.Bchain,
            row.antigen_chain,
            row.mhc_chain1,
            row.mhc_chain2,
            row.mhc_type,
            row.imgt_file_path,
            row.raw_file_path,
        ),
        axis=1,
    )

    return (
        df[(structure_status != 'missing') & (structure_status != 'unknown')].copy(),
        structure_status[(structure_status != 'missing') & (structure_status != 'unknown')].copy(),
    )


def add_cdr_sequences(pdb_id: str, alpha_chain_id: str, beta_chain_id: str, imgt_file_path: str) -> pd.Series:
    """Add CDR sequences from structures."""
    structure = PDBParser().get_structure(pdb_id, imgt_file_path)

    sequences = {
        f'cdr_{chain_type[0]}{cdr_num}': get_sequence(structure, chain_id, numbering)
        for chain_type, chain_id in (('alpha_chain', alpha_chain_id), ('beta_chain', beta_chain_id))
        for cdr_num, numbering in ((1, IMGT_CDR1), (2, IMGT_CDR2), (3, IMGT_CDR3))
    }

    return pd.Series(sequences).sort_index()


def add_peptide_sequences(pdb_id: str, antigen_chain_id: str, imgt_file_path: str) -> str:
    """Add peptide sequences from structures."""
    structure = PDBParser().get_structure(pdb_id, imgt_file_path)

    return get_sequence(structure, antigen_chain_id)


def add_mhc_tcr_contact_pseudo_sequences(
    pdb_id: str,
    mhc_type: str,
    mhc_chain_1_id: str,
    mhc_chain_2_id: str,
    stcrdab_path: str,
    mhc_tcr_contact_residues: set[int] | tuple[set[int], set[int]],
) -> str:
    """Add MHC-TCR contact pseudo sequences."""
    structure = PDBParser().get_structure(pdb_id, os.path.join(stcrdab_path, 'imgt', pdb_id + '.pdb'))

    if mhc_type == 'MH1':
        return get_sequence(structure, mhc_chain_1_id, mhc_tcr_contact_residues)

    if mhc_type == 'MH2':
        return get_sequence(structure, mhc_chain_1_id, mhc_tcr_contact_residues[0]) + get_sequence(
            structure, mhc_chain_2_id, mhc_tcr_contact_residues[1]
        )

    msg = f'Invalid MHC type: {mhc_type}. Type must be MH1 or MH2.'
    raise ValueError(msg)


def remove_similar_structures(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Remove structures with the same CDR and peptide sequences within the RMSD threshold.

    The highest resolution structure will be kept.
    """
    output_dfs = []

    for (cdr_sequence, peptide_sequence, mhc_type), group in df.groupby(
        ['collated_cdrs', 'peptide_sequence', 'mhc_type']
    ):
        if len(group) == 1:
            output_dfs.append(group)
            continue

        logger.debug('Screening TCR: %s, peptide: %s, MHC: %s', cdr_sequence, peptide_sequence, mhc_type)

        sorted_group = group.sort_values(
            ['resolution', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2'],
        )

        pdb_parser = PDBParser(QUIET=True)

        structures = []
        chain_maps = []

        for _, row in sorted_group.iterrows():
            chain_map = {'alpha_chain': row.Achain, 'beta_chain': row.Bchain, 'antigen_chain': row.antigen_chain}

            if mhc_type == 'MH1':
                chain_map['mhc_chain1'] = row.mhc_chain1

            elif mhc_type == 'MH2':
                chain_map['mhc_chain1'] = row.mhc_chain1
                chain_map['mhc_chain2'] = row.mhc_chain2

            else:
                msg = f"Invalid MHC type {mhc_type}. MHC type must be 'MH1' or 'MH2'"
                raise ValueError(msg)

            logger.debug('Collecting PDB ID: %s and extracting chains %s', row.pdb, '-'.join(chain_map.values()))

            structure = pdb_parser.get_structure('', row.imgt_file_path)
            structure = extract_chains(structure, chain_map.values())

            structures.append(structure)
            chain_maps.append(chain_map)

        distance_matrix = compute_structural_distances(structures, chain_maps, mhc_type)
        clusters = (
            AgglomerativeClustering(
                metric='precomputed', distance_threshold=threshold, linkage='single', n_clusters=None
            )
            .fit(distance_matrix)
            .labels_
        )
        clusters = pd.Series(clusters, index=sorted_group.index)
        output_dfs.append(sorted_group[~clusters.duplicated()])

    return pd.concat(output_dfs)


def split_groups_multiple_proportions(idx_sizes, proportions, exclude_idxs=None, exclude_split=None):
    total_sum = sum([size for _, size in idx_sizes])
    target_sums = [total_sum * prop for prop in proportions]

    random.shuffle(idx_sizes)

    groups = [[] for _ in range(len(proportions))]
    sum_groups = [0] * len(proportions)

    if exclude_idxs is not None:  # add all excluded IDs to selected split
        for idx, size in idx_sizes:
            if idx in exclude_idxs and exclude_split is not None:
                groups[exclude_split].append(idx)
                sum_groups[exclude_split] += size

    for idx, size in idx_sizes:
        if exclude_idxs and idx in exclude_idxs:
            continue

        # Find the group with the smallest current sum and add the number
        min_sum_index = min(range(len(proportions)), key=lambda i: sum_groups[i])

        if sum_groups[min_sum_index] + size <= target_sums[min_sum_index]:
            groups[min_sum_index].append(idx)
            sum_groups[min_sum_index] += size

        else:
            # If adding the number exceeds the target sum, add to the next group
            for i in range(len(proportions)):
                if i != min_sum_index and sum_groups[i] + size <= target_sums[i]:
                    groups[i].append(idx)
                    sum_groups[i] += size
                    break
            else:  # If it doesn't fit anywhere, add it to the original trial
                groups[min_sum_index].append(idx)
                sum_groups[min_sum_index] += size

    return groups


def merge_groups(groups):
    merged_groups = []

    for group in groups:
        merged = False

        for i, merged_group in enumerate(merged_groups):
            if len(group & merged_group) > 0:
                merged_groups[i] = group | merged_group
                merged = True
                break

        if not merged:
            merged_groups.append(group)

    return merged_groups


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    if args.seed:
        random.seed(args.seed)

    if args.fix_structures_missing_residues:
        fix_dir = tempfile.TemporaryDirectory()
        fix_dir_name = fix_dir.name

    else:
        fix_dir_name = None

    stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab, 'db_summary.dat'), delimiter='\t')
    stcrdab_summary['imgt_file_path'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'imgt', pdb_id + '.pdb'),
    )
    stcrdab_summary['raw_file_path'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'raw', pdb_id + '.pdb'),
    )

    logger.info('Selecting structures...')
    selected_structures = stcrdab_summary
    selected_structures = selected_structures.query('TCRtype in @args.tcr_types')
    selected_structures = selected_structures.query('mhc_type in @args.mhc_types')
    selected_structures = selected_structures.query('antigen_type in @args.antigen_types')
    selected_structures = selected_structures.copy()

    logger.info('Screening Quality...')
    selected_structures['resolution'] = pd.to_numeric(selected_structures['resolution'], errors='coerce')
    selected_structures = selected_structures.query('resolution <= @args.resolution_cutoff')

    if args.remove_structures_missing_residues:
        logger.debug('Removing structures missing residues')
        selected_structures, structure_status = screen_for_missing_residues(
            selected_structures,
            fix_structures=args.fix_structures_missing_residues,
            fix_dir=fix_dir_name,
        )

        if args.fix_structures_missing_residues:
            selected_structures.loc[structure_status == 'fixed', 'imgt_file_path'] = selected_structures[
                structure_status == 'fixed'
            ]['imgt_file_path'].map(
                lambda path: os.path.join(fix_dir_name, os.path.basename(path)),
            )

    logger.info('Getting sequence information...')
    cdr_sequences = selected_structures.apply(
        lambda row: add_cdr_sequences(row.pdb, row.Achain, row.Bchain, row.imgt_file_path),
        axis=1,
    )

    peptide_sequences = selected_structures.apply(
        lambda row: add_peptide_sequences(row.pdb, row.antigen_chain, row.imgt_file_path),
        axis=1,
    )
    peptide_sequences.name = 'peptide_sequence'

    mhc_tcr_contacts_available = args.mhc_class_I_tcr_contact_residues or (
        args.mhc_class_II_alpha_chain_tcr_contact_residues and args.mhc_class_II_beta_chain_tcr_contact_residues
    )
    if mhc_tcr_contacts_available:
        mhc_i_tcr_contact_residues_range = {
            int(''.join([character for character in seq_id if character.isnumeric()]))
            for seq_id in args.mhc_class_I_tcr_contact_residues
        }

        mhc_ii_tcr_contact_residues_range = (
            {
                int(''.join([character for character in seq_id if character.isnumeric()]))
                for seq_id in args.mhc_class_II_alpha_chain_tcr_contact_residues
            },
            {
                int(''.join([character for character in seq_id if character.isnumeric()]))
                for seq_id in args.mhc_class_II_beta_chain_tcr_contact_residues
            },
        )

        mhc_tcr_contact_pseudo_sequences = selected_structures.apply(
            lambda row: add_mhc_tcr_contact_pseudo_sequences(
                row.pdb,
                row.mhc_type,
                row.mhc_chain1,
                row.mhc_chain2,
                args.stcrdab,
                mhc_i_tcr_contact_residues_range if row.mhc_type == 'MH1' else mhc_ii_tcr_contact_residues_range,
            ),
            axis=1,
        )
        mhc_tcr_contact_pseudo_sequences.name = 'mhc_tcr_contact_pseudo_sequence'

        selected_structures = pd.concat(
            [
                selected_structures,
                cdr_sequences,
                peptide_sequences.to_frame(),
                mhc_tcr_contact_pseudo_sequences.to_frame(),
            ],
            axis='columns',
        )

    else:
        selected_structures = pd.concat(
            [selected_structures, cdr_sequences, peptide_sequences.to_frame()], axis='columns'
        )

    selected_structures['collated_cdrs'] = (
        selected_structures['cdr_a1']
        + '-'
        + selected_structures['cdr_a2']
        + '-'
        + selected_structures['cdr_a3']
        + '-'
        + selected_structures['cdr_b1']
        + '-'
        + selected_structures['cdr_b2']
        + '-'
        + selected_structures['cdr_b3']
    )

    if args.structural_similarity_cutoff:
        logger.info('Removing structures within %.2f Å RMSD', args.structural_similarity_cutoff)
        selected_structures = remove_similar_structures(selected_structures, args.structural_similarity_cutoff)

    logger.info(
        'Splitting data accoding to partions (Train: %.2f, Validation %.2f, and Test %.2f)',
        args.train_split,
        args.validation_split,
        args.test_split,
    )

    peptide_groups = selected_structures.groupby('peptide_sequence')

    logger.debug('Merging peptide groups with common TCRs.')
    merge_matrix = np.zeros((len(peptide_groups), len(peptide_groups)))
    for i, (_, group_i) in enumerate(peptide_groups):
        for j, (_, group_j) in enumerate(peptide_groups):
            group_i_seqs = set(group_i['collated_cdrs'].tolist())
            group_j_seqs = set(group_j['collated_cdrs'].tolist())

            if len(group_i_seqs & group_j_seqs) > 0:
                merge_matrix[i, j] = 1

    group_indicies = np.arange(len(peptide_groups))
    groups = [set(group_indicies[row > 0]) for row in merge_matrix]

    separated_groups = merge_groups(groups)

    separated_groups_data = [
        pd.concat([list(peptide_groups)[idx][1] for idx in group], axis=0) for group in separated_groups
    ]

    logger.debug('Adding group IDs')
    for group_id, group in enumerate(separated_groups_data, 1):
        group['group_id'] = group_id

    exclude_group_idxs = (
        [idx for idx, group in enumerate(separated_groups_data) if group['pdb'].isin(args.pdb_ids_to_exclude).any()]
        if args.pdb_ids_to_exclude is not None
        else None
    )

    separated_idx_size = [(idx, len(df)) for idx, df in enumerate(separated_groups_data)]

    train_idxs, val_idxs, test_idxs = split_groups_multiple_proportions(
        separated_idx_size,
        (args.train_split, args.validation_split, args.test_split),
        exclude_group_idxs,
        exclude_split=0,
    )  # training split

    dataset = pd.DataFrame()
    for split_name, idxs in (('train', train_idxs), ('validation', val_idxs), ('test', test_idxs)):
        if len(idxs) == 0:
            logger.warning('No data in %s split', split_name)
            continue

        split_data = pd.concat([df for idx, df in enumerate(separated_groups_data) if idx in idxs], axis=0)
        split_data['split'] = split_name

        dataset = pd.concat([dataset, split_data])

    dataset = dataset.fillna('')
    dataset = dataset.reset_index()

    logger.info('Outputing structures...')
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    output_columns = [
        'pdb',
        'Achain',
        'Bchain',
        'antigen_chain',
        'mhc_chain1',
        'mhc_chain2',
        'mhc_type',
        'peptide_sequence',
        'collated_cdrs',
    ]

    if mhc_tcr_contacts_available:
        output_columns.append('mhc_tcr_contact_pseudo_sequence')

    output_columns += ['group_id', 'split']

    dataset[output_columns].to_csv(os.path.join(args.output, 'stcrdab_split.csv'), index=False)

    pdb_parser = PDBParser()
    for _, row in dataset.iterrows():
        structure = pdb_parser.get_structure(row.pdb, row.imgt_file_path)
        output_name = f'{row.pdb}_{row.Achain}{row.Bchain}{row.antigen_chain}{row.mhc_chain1}{row.mhc_chain2}.pdb'
        output_chains = [row.Achain, row.Bchain, row.antigen_chain, row.mhc_chain1, row.mhc_chain2]

        if args.crop_structures:
            logger.debug('Cropping structure')
            match row.mhc_type:
                case 'MH1':
                    structure = crop_structure(structure, (row.Achain, row.Bchain), (row.mhc_chain1,), row.mhc_type)

                    output_chains = [row.Achain, row.Bchain, row.antigen_chain, row.mhc_chain1]

                case 'MH2':
                    structure = crop_structure(
                        structure, (row.Achain, row.Bchain), (row.mhc_chain1, row.mhc_chain2), row.mhc_type
                    )

        io = PDBIO()
        io.set_structure(structure)
        io.save(os.path.join(args.output, output_name), SelectChains(*output_chains))

    if args.fix_structures_missing_residues:
        fix_dir.cleanup()


if __name__ == '__main__':
    main()
