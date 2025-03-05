"""Process contact information to use in downstream machine learning tasks."""

import argparse
import json
import logging
import sys

import h5py
import numpy as np
import pandas as pd

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import MHC_I_IMGT_BETA_HELIX_START, assign_helix

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('contacts', help='path to input contact data file')
parser.add_argument(
    '--mhc-pseudo-sequence-imgt-numbers',
    required=True,
    help='path to mhc pseudo sequence imgt numbers (JSON format)',
)
parser.add_argument('--output', '-o', required=True, help='Path to output h5 file')

data_group = parser.add_argument_group('Data')
data_group.add_argument(
    '--mhc-types',
    nargs='+',
    choices=['MH1', 'MH2'],
    default=['MH1', 'MH2'],
    help="TCR:pMHC types to include (default: ['MH1', 'MH2'])",
)
data_group.add_argument(
    '--smooth',
    action='store_true',
    help='Use gaussian to smooth contact maps between TCRs and peptides',
)
data_group.add_argument('--pad', action='store_true', help='Pad maps using the provided lengths')
data_group.add_argument(
    '--cdr1-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 1a sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr2-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 2a sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr3-alpha-length',
    default=24,
    type=int,
    help='Length to pad CDR 3a sequences to (Default: 24)',
)
data_group.add_argument(
    '--cdr1-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 1b sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr2-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 2b sequences to (Default: 8)',
)
data_group.add_argument(
    '--cdr3-beta-length',
    default=24,
    type=int,
    help='Length to pad CDR 3b sequences to (Default: 24)',
)
data_group.add_argument(
    '--peptide-length',
    default=12,
    type=int,
    help='Length to pad peptide sequences to (Default: 12)',
)


add_logging_arguments(parser)


def gaussian_2d(values: np.ndarray, means: np.ndarray, cov: np.ndarray) -> np.ndarray:
    return (1 / np.sqrt(((2 * np.pi) ** 2) * np.linalg.det(cov))) * np.exp(
        -1 / 2 * (values - means).T @ np.linalg.inv(cov) @ (values - means)
    )


def split_resi(resi: str) -> tuple[int, str]:
    return (
        int(''.join([char for char in resi if char.isnumeric()])),
        '' if resi[-1].isnumeric() else resi[-1],
    )


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading contact data')
    contacts = pd.read_csv(args.contacts)

    logger.debug('Formatting TCR data')
    contacts['chain_type_tcr'] = contacts['cdr_type'].str[4:]
    contacts['cdr'] = contacts['cdr_type'].str[3].apply(int)

    logger.debug('Formatting MHC data')
    with open(args.mhc_pseudo_sequence_imgt_numbers, 'r') as fh:
        mhc_pseudo_sequences = json.load(fh)

    mhc_pseudo_pos = [
        (helix[0], *split_resi(resi)) for helix in ('alpha', 'beta') for resi in mhc_pseudo_sequences[helix]
    ]

    contacts['mhc_pos'] = None
    contacts.loc[contacts['resi_pmhc'].notna(), 'mhc_pos'] = contacts.loc[
        contacts['resi_pmhc'].notna(), 'resi_pmhc'
    ].map(split_resi)

    contacts['residue_seq_id_pmhc'] = None
    contacts.loc[contacts['resi_pmhc'].notna(), 'residue_seq_id_pmhc'] = contacts[contacts['resi_pmhc'].notna()][
        'mhc_pos'
    ].map(lambda resi: resi[0])

    contacts['imgt_pmhc'] = None
    shortened_resis = contacts[contacts['resi_pmhc'].notna()]['mhc_pos'].map(
        lambda resi: resi[0] - MHC_I_IMGT_BETA_HELIX_START if resi[0] > MHC_I_IMGT_BETA_HELIX_START else resi[0],
    )
    insert_codes = contacts[contacts['resi_pmhc'].notna()]['mhc_pos'].map(lambda resi: resi[1])
    mhc_helices = contacts[contacts['resi_pmhc'].notna()].apply(
        lambda row: assign_helix(row.mhc_type, row.chain_type_pmhc, row.resi_pmhc)[0],
        axis=1,
    )
    contacts.loc[contacts['resi_pmhc'].notna(), 'mhc_pos'] = pd.concat(
        [mhc_helices.to_frame(), shortened_resis.to_frame(), insert_codes.to_frame()],
        axis=1,
    ).apply(tuple, axis=1)

    if args.pad:
        contacts['cdr_pos'] = None

        for chain_type, cdr_num, cdr_length in (
            ('alpha', 1, args.cdr1_alpha_length),
            ('alpha', 2, args.cdr2_alpha_length),
            ('alpha', 3, args.cdr3_alpha_length),
            ('beta', 1, args.cdr1_beta_length),
            ('beta', 2, args.cdr2_beta_length),
            ('beta', 3, args.cdr3_beta_length),
        ):
            selection = (
                (contacts['chain_type_tcr'] == chain_type)
                & (contacts['cdr'] == cdr_num)
                & contacts['relative_pos_centre_tcr'].notna()
            )
            contacts.loc[selection, 'cdr_pos'] = contacts.loc[selection, 'relative_pos_centre_tcr'].map(
                lambda pos, cdr_length=cdr_length: pos + cdr_length // 2 - 1,
            )

        contacts['peptide_pos'] = None
        contacts.loc[contacts['relative_pos_centre_pmhc'].notna(), 'peptide_pos'] = contacts[
            contacts['relative_pos_centre_pmhc'].notna()
        ]['relative_pos_centre_pmhc'].map(
            lambda pos: pos + args.peptide_length // 2 - 1,
        )

    else:
        contacts['cdr_pos'] = contacts['relative_pos_centre_tcr']
        contacts['peptide_pos'] = contacts['relative_pos_centre_pmhc']

    contact_maps = {}

    logger.info('Compiling TCR:peptide contact maps')
    contact_maps['peptide'] = {}
    tcr_peptide_data = contacts[
        (contacts['mhc_type'].isin(args.mhc_types)) & (contacts['chain_type_pmhc'] == 'antigen_chain')
    ]

    y_labels = sorted(tcr_peptide_data['cdr_pos'].unique())
    x_labels = sorted(tcr_peptide_data['peptide_pos'].unique())

    for chain_type, cdr_num, cdr_length in (
        ('alpha', 1, args.cdr1_alpha_length),
        ('alpha', 2, args.cdr2_alpha_length),
        ('alpha', 3, args.cdr3_alpha_length),
        ('beta', 1, args.cdr1_beta_length),
        ('beta', 2, args.cdr2_beta_length),
        ('beta', 3, args.cdr3_beta_length),
    ):
        cdr_type = f'cdr{cdr_num}_{chain_type}'
        logger.debug('Compiling %s:peptide contact map', cdr_type)

        cdr_data = tcr_peptide_data[
            (tcr_peptide_data['chain_type_tcr'] == chain_type) & (tcr_peptide_data['cdr'] == cdr_num)
        ]

        if args.pad:
            logger.debug('Padding and cropping contact map')
            y_labels = np.arange(cdr_length)
            x_labels = np.arange(args.peptide_length)

            cdr_data = cdr_data[
                (cdr_data['cdr_pos'] >= 0)
                & (cdr_data['cdr_pos'] < cdr_length)
                & (cdr_data['peptide_pos'] >= 0)
                & (cdr_data['peptide_pos'] < args.peptide_length)
            ]

        if args.smooth:
            logger.debug('Smoothing map with a gaussian distribution')
            group_max = cdr_data['proportion'].max()

            x_dist = cdr_data.groupby(['peptide_pos'])['proportion'].sum().reset_index()
            x_mu = np.average(x_dist['peptide_pos'], weights=x_dist['proportion'])
            x_var = np.average((x_dist['peptide_pos'] - x_mu) ** 2, weights=x_dist['proportion'])

            y_dist = cdr_data.groupby(['cdr_pos'])['proportion'].sum().reset_index()
            y_mu = np.average(y_dist['cdr_pos'], weights=y_dist['proportion'])
            y_var = np.average((y_dist['cdr_pos'] - y_mu) ** 2, weights=y_dist['proportion'])

            means = np.array([x_mu, y_mu])
            cov = np.array([[x_var, 0], [0, y_var]])

            gaussian = []
            for x in x_labels:
                for y in y_labels:
                    gauss_xy = gaussian_2d(np.array([x, y]), means, cov)
                    gaussian.append([x, y, gauss_xy])

            gaussian = pd.DataFrame(gaussian, columns=['peptide_pos', 'cdr_pos', 'gaussian'])
            gaussian['gaussian'] = gaussian['gaussian'] * (group_max / gaussian['gaussian'].max())
            contact_map = gaussian.pivot_table(values='gaussian', index='cdr_pos', columns='peptide_pos')

        else:
            contact_map = cdr_data.pivot_table(index='cdr_pos', columns='peptide_pos', values='proportion')

            logger.debug('Filling missing values in contact map')
            contact_map[[label for label in x_labels if label not in contact_map.columns]] = np.nan
            missing_values = pd.DataFrame(
                np.nan,
                index=pd.Series(
                    [label for label in y_labels if label not in contact_map.index],
                    name=contact_map.index.name,
                ),
                columns=contact_map.columns,
            )
            if not missing_values.empty:
                contact_map = pd.concat([contact_map, missing_values])

            contact_map = contact_map.fillna(0)

        contact_map = contact_map[sorted(contact_map.columns)].sort_index(ascending=False)

        contact_maps['peptide'][cdr_type] = contact_map.to_numpy()

    logger.info('Normalising TCR:peptide contact maps')
    total = sum(contact_map.sum() for contact_map in contact_maps['peptide'].values())
    contact_maps['peptide'] = {
        cdr_type: contact_map / total for cdr_type, contact_map in contact_maps['peptide'].items()
    }

    logger.info('Compiling TCR:MHC contact maps')
    contact_maps['mhc_pseudo'] = {}
    tcr_mhc_data = contacts[
        contacts['mhc_type'].isin(args.mhc_types) & (contacts['chain_type_pmhc'] != 'antigen_chain')
    ]

    y_labels = sorted(tcr_mhc_data['cdr_pos'].unique())
    x_labels = sorted(tcr_mhc_data['mhc_pos'].unique())

    for chain_type, cdr_num, cdr_length in (
        ('alpha', 1, args.cdr1_alpha_length),
        ('alpha', 2, args.cdr2_alpha_length),
        ('alpha', 3, args.cdr3_alpha_length),
        ('beta', 1, args.cdr1_beta_length),
        ('beta', 2, args.cdr2_beta_length),
        ('beta', 3, args.cdr3_beta_length),
    ):
        cdr_type = f'cdr{cdr_num}_{chain_type}'
        logger.debug('Compiling %s:MHC contact map', cdr_type)

        cdr_data = tcr_mhc_data[(tcr_mhc_data['chain_type_tcr'] == chain_type) & (tcr_mhc_data['cdr'] == cdr_num)]

        if args.pad:
            logger.debug('Padding and cropping contact map')
            y_labels = np.arange(cdr_length)
            x_labels = mhc_pseudo_pos

            cdr_data = cdr_data[
                (cdr_data['cdr_pos'] >= 0)
                & (cdr_data['cdr_pos'] < cdr_length)
                & cdr_data['mhc_pos'].isin(mhc_pseudo_pos)
            ]

        contact_map = cdr_data.pivot_table(index='cdr_pos', columns='mhc_pos', values='proportion')

        logger.debug('Filling missing values in contact map')
        contact_map[[label for label in x_labels if label not in contact_map.columns]] = np.nan
        missing_values = pd.DataFrame(
            np.nan,
            index=pd.Series(
                [label for label in y_labels if label not in contact_map.index],
                name=contact_map.index.name,
            ),
            columns=contact_map.columns,
        )
        if not missing_values.empty:
            contact_map = pd.concat([contact_map, missing_values])

        contact_map = contact_map.fillna(0)

        contact_map = contact_map[sorted(contact_map.columns)].sort_index(ascending=False)

        contact_maps['mhc_pseudo'][cdr_type] = contact_map.to_numpy()

    logger.info('Normalising TCR:MHC contact maps')
    total = sum(contact_map.sum() for contact_map in contact_maps['mhc_pseudo'].values())
    contact_maps['mhc_pseudo'] = {
        cdr_type: contact_map / total for cdr_type, contact_map in contact_maps['mhc_pseudo'].items()
    }

    logger.info('Outputting contact maps to %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        for group_name, cdr_data in contact_maps.items():
            group = fh.create_group(group_name)

            for cdr_type, values in cdr_data.items():
                group[cdr_type] = values


if __name__ == '__main__':
    main()
