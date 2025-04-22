"""Count the types of interacting resiudes from the PPI3D dataset https://bioinformatics.lt/ppi3d/."""

import argparse
import logging
import sys
import tempfile

import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData
from scipy.spatial import KDTree

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.structure import bio_to_pandas

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input', help='Path to the csv with entry information from PPI3D.')
parser.add_argument('--output', '-o', required=True, help='Path to output the csv with contacting residue information.')
parser.add_argument('--seed', default=None, type=int, help='Random seed for data splitting (Default: None).')
parser.add_argument(
    '--sample-size',
    default=None,
    type=int,
    help=(
        'Sample size for random samples (Default: None). If this parameter is left as None, the whole input dataset '
        'will be used.'
    ),
)
parser.add_argument(
    '--normalise',
    action='store_true',
    help='Add this flag to normalise the final output over the total number of interacting residues.',
)
parser.add_argument(
    '--contact-distance',
    type=float,
    default=5.0,
    help='Maximum distance to be considered a contact (Default: 5 Å).',
)
parser.add_argument('--separate-interaction-types', action='store_true', help='')

add_logging_arguments(parser)


def find_contacting_residues(
    struct_df1: pd.DataFrame, struct_df2: pd.DataFrame, contact_distance: float
) -> pd.DataFrame:
    """Find contacting residues between two structures in pandas format."""
    coords1 = struct_df1[['pos_x', 'pos_y', 'pos_z']].to_numpy()
    coords2 = struct_df2[['pos_x', 'pos_y', 'pos_z']].to_numpy()

    tree1 = KDTree(coords1)
    tree2 = KDTree(coords2)

    pairs = tree1.query_ball_tree(tree2, contact_distance)

    contact_rows = []
    for i, indices in enumerate(pairs):
        for j in indices:
            row_1 = struct_df1.iloc[i]
            row_2 = struct_df2.iloc[j]
            contact_rows.append(
                {
                    **{f'{col}_1': row_1[col] for col in struct_df1.columns},
                    **{f'{col}_2': row_2[col] for col in struct_df2.columns},
                }
            )

    contacts = pd.DataFrame(contact_rows)

    return contacts


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    if args.seed:
        logger.info('Setting seed to %d', args.seed)
        rng = np.random.default_rng(args.seed)

    else:
        rng = np.random.default_rng()

    logger.info('Loading metadata')
    ppi_3d_summary = pd.read_csv(args.input)

    if args.separate_interaction_types:
        interaction_types = sorted(ppi_3d_summary.filter(regex='_interaction$').columns.tolist())
        ppi_3d_summary[['interaction_type']] = pd.from_dummies(ppi_3d_summary[interaction_types])

    if args.sample_size:
        logger.info('Sampling %d samples from the dataset', args.sample_size)
        ppi_3d_summary = ppi_3d_summary.sample(args.sample_size, random_state=rng).reset_index()

    num_ppis = len(ppi_3d_summary)

    logger.info('Counting residue interactions')
    amino_acid_names = sorted([name.upper() for name in IUPACData.protein_letters_3to1])

    if args.separate_interaction_types:
        residue_counts = pd.Series(
            0,
            index=pd.MultiIndex.from_product(
                [interaction_types, amino_acid_names, amino_acid_names],
                names=['interaction_type', 'residue_name_1', 'residue_name_2'],
            ),
            name='count',
        )

    else:
        residue_counts = pd.Series(
            0,
            index=pd.MultiIndex.from_product(
                [amino_acid_names, amino_acid_names], names=['residue_name_1', 'residue_name_2']
            ),
            name='count',
        )

    for i, (_, ppi_entry) in enumerate(ppi_3d_summary.iterrows(), 1):
        logger.debug('Working on %s - %d of %d', ppi_entry.pdb_id, i, num_ppis)

        logger.debug('Downloading data')
        req = requests.get(ppi_entry.download_url, timeout=60)

        with tempfile.NamedTemporaryFile('w+') as fh:
            fh.write(req.text)

            pdb_parser = PDBParser(QUIET=True)
            structure = pdb_parser.get_structure('', fh.name)

        logger.debug('Cleaning structure')
        structure_df = bio_to_pandas(structure)
        structure_df_clean = structure_df[structure_df['record_type'] == 'ATOM']
        structure_df_clean = structure_df_clean[structure_df_clean['residue_name'].isin(amino_acid_names)]

        subunit_1_df = structure_df_clean[structure_df_clean['chain_id'] == 'A']
        subunit_2_df = structure_df_clean[structure_df_clean['chain_id'] == 'B']

        logger.debug('Computing interface residues at a distance of <%d Å', args.contact_distance)
        contacts = find_contacting_residues(subunit_1_df, subunit_2_df, args.contact_distance)
        contacts = contacts.drop_duplicates(
            [
                'chain_id_1',
                'residue_seq_id_1',
                'residue_insert_code_1',
                'chain_id_2',
                'residue_seq_id_2',
                'residue_insert_code_2',
            ]
        )
        logger.debug('%d residue contact(s) found', len(contacts))

        logger.debug('Updating contact counts')
        if args.separate_interaction_types:
            residue_counts[ppi_entry['interaction_type']] = (
                residue_counts[ppi_entry['interaction_type']]
                .add(contacts.value_counts(['residue_name_1', 'residue_name_2']), fill_value=0)
                .astype(int)
            )

        else:
            residue_counts = residue_counts.add(
                contacts.value_counts(['residue_name_1', 'residue_name_2']), fill_value=0
            ).astype(int)

    if args.normalise:
        logger.info('Normalising counts')
        residue_counts = residue_counts / residue_counts.sum()
        residue_counts.name = 'proportion'

    logger.info('Outputting residue counts')
    residue_counts.to_frame().to_csv(args.output)


if __name__ == '__main__':
    main()
