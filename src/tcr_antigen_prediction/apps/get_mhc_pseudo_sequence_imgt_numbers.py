"""Get the relevant IMGT numbers for to create MHC Pseudo sequences from TCR contact maps."""

import argparse
import collections
import json
import logging
import re
import sys

import pandas as pd

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('contact_map', help='Path to the contact map csv.')
parser.add_argument('--output', '-o', help='Path to the output json file.')
parser.add_argument(
    '--probability-filter',
    type=float,
    default=0.01,
    help='Cutoff to keep positions above this proportion of contacts (Default: 0.01).',
)

add_logging_arguments(parser)

MHC_I_IMGT_BETA_HELIX_START = 1000


def assign_helix(mhc_type: str, chain_type: str, resi: str) -> str:
    """Assign an MHC residue as being either part of the 'alpha' helix or 'beta' helix.

    The assignment happens regardless of class I versus class II.

    """
    match mhc_type:
        case 'MH1':
            return (
                'alpha'
                if int(''.join([char for char in resi if char.isnumeric()])) < MHC_I_IMGT_BETA_HELIX_START
                else 'beta'
            )

        case 'MH2':
            return 'alpha' if chain_type == 'mhc_chain1' else 'beta'


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading data')
    tcr_pmhc_contacts = pd.read_csv(args.contact_map)
    tcr_mhc_contacts = tcr_pmhc_contacts.query("chain_type_pmhc != 'antigen_chain'").copy()

    logger.info('Annotating data')
    tcr_mhc_contacts['helix'] = tcr_mhc_contacts.apply(
        lambda row: assign_helix(row.mhc_type, row.chain_type_pmhc, row.resi_pmhc),
        axis=1,
    )
    mhc_1_chain_2_pattern = re.compile(r'10\d\d[A-Z]?')
    tcr_mhc_contacts['resi_pmhc_imgt'] = tcr_mhc_contacts['resi_pmhc'].map(
        lambda resi: resi[2:] if re.match(mhc_1_chain_2_pattern, resi) else resi,
    )

    logger.info('Re-normalising MHC data')
    mhc_contacts = (
        tcr_mhc_contacts.groupby(
            ['mhc_type', 'chain_type_pmhc', 'helix', 'resi_pmhc', 'resi_pmhc_imgt'],
        )['proportion']
        .sum()
        .reset_index()
    )

    mhc_contacts['proportion'] = mhc_contacts.groupby('mhc_type')['proportion'].transform(
        lambda proportion: proportion / proportion.sum(),
    )

    mhc_contacts = mhc_contacts.groupby(['helix', 'resi_pmhc_imgt'])['proportion'].sum().reset_index()
    mhc_contacts['proportion'] /= mhc_contacts['proportion'].sum()

    logger.info('Filtering contacts below probability filter (%.2f)', args.probability_filter)
    mhc_contacts_filtered = mhc_contacts[mhc_contacts['proportion'] > args.probability_filter]

    logger.info('Sorting data by chain and IMGT number')
    mhc_contacts_filtered_sorted = mhc_contacts_filtered.copy()
    mhc_contacts_filtered_sorted['residue_seq_id'] = mhc_contacts_filtered_sorted['resi_pmhc_imgt'].map(
        lambda resi: int(''.join(char for char in resi if char.isnumeric())),
    )
    mhc_contacts_filtered_sorted['residue_insert_code'] = mhc_contacts_filtered_sorted['resi_pmhc_imgt'].map(
        lambda resi: ''.join(char for char in resi if not char.isnumeric()),
    )

    mhc_contacts_filtered_sorted = mhc_contacts_filtered_sorted.sort_values(
        ['helix', 'residue_seq_id', 'residue_insert_code'],
    )
    mhc_contacts_filtered_sorted = mhc_contacts_filtered_sorted.drop(columns=['residue_seq_id', 'residue_insert_code'])

    logger.info('Collecting IMGT residues')
    mhc_pseudo_seq_positions = collections.defaultdict(list)
    for _, row in mhc_contacts_filtered_sorted[['helix', 'resi_pmhc_imgt']].iterrows():
        mhc_pseudo_seq_positions[row['helix']].append(row['resi_pmhc_imgt'])

    logger.info('Outputting data')
    with open(args.output, 'w') as fh:
        json.dump(mhc_pseudo_seq_positions, fh)


if __name__ == '__main__':
    main()
