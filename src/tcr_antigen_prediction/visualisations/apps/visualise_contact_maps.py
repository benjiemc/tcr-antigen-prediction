"""Visualise TCR:pMHC Contact Maps.

CDR loops and peptides are standardised by position relative to the apex of the loop. In cases where the entities length
is even, the position to the left of the middle of the loop is taken as the centre.

"""

import argparse
import json
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import MHC_I_IMGT_BETA_HELIX_START, assign_helix

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('contacts', help='path to csv with contact data')
parser.add_argument('--output', '-o', help='path to output plot')
parser.add_argument('--smooth', action='store_true', help='fit data to a 2D gaussian to smooth')
parser.add_argument('--normalise', action='store_true', help='normalise data for all proportions to sum to 1')
parser.add_argument('--separate-plots', action='store_true', help='save a separate plot for each CDR')

padding_group = parser.add_argument_group('Padding')
padding_group.add_argument('--pad', action='store_true', help='centre pad sequences to specified length')
padding_group.add_argument('--cdr-1-length', type=int, default=8, help='pad length for CDR1s (Default: 8)')
padding_group.add_argument('--cdr-2-length', type=int, default=8, help='pad length for CDR2s (Default: 8)')
padding_group.add_argument('--cdr-3-length', type=int, default=24, help='pad length for CDR3s (Default: 24)')
padding_group.add_argument('--peptide-length', type=int, default=12, help='pad length for peptides (Default: 12)')
padding_group.add_argument('--mhc-pseudo-seq-positions', help='MHC pseudo sequence positions to use for padding MHCs')

selection_group = parser.add_argument_group('Selection')
selection_group.add_argument(
    '--interaction',
    choices=['peptide', 'MHC'],
    default='peptide',
    help="choice of entity to plot interactions. eg: TCR:peptide or TCR:MHC (Default: 'peptide')",
)
selection_group.add_argument(
    '--mhc-types',
    nargs='+',
    default=['MH1', 'MH2'],
    help="MHC types to include (Default: 'MH1' and 'MH2')",
)

add_logging_arguments(parser)


def gaussian_2d(values: np.ndarray, means: np.ndarray, cov: np.ndarray) -> np.ndarray:
    return (1 / np.sqrt(((2 * np.pi) ** 2) * np.linalg.det(cov))) * np.exp(
        -1 / 2 * (values - means).T @ np.linalg.inv(cov) @ (values - means)
    )


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading contacts from %s', args.contacts)
    contact_map = pd.read_csv(args.contacts)

    contact_map['chain_type_tcr'] = contact_map['cdr_type'].str[4:]
    contact_map['cdr'] = contact_map['cdr_type'].str[3].apply(int)

    contact_map['residue_seq_id_pmhc'] = None
    contact_map.loc[contact_map['resi_pmhc'].notna(), 'residue_seq_id_pmhc'] = contact_map[
        contact_map['resi_pmhc'].notna()
    ]['resi_pmhc'].map(
        lambda resi: int(''.join([char for char in resi if char.isnumeric()])),
    )

    if args.interaction == 'MHC':
        contact_map.loc[pd.notna(contact_map['resi_pmhc']), 'resi_pmhc'] = contact_map.loc[
            pd.notna(contact_map['resi_pmhc']), 'resi_pmhc'
        ].map(
            lambda resi: (
                int(''.join([char for char in resi if char.isnumeric()])),
                '' if resi[-1].isnumeric() else resi[-1],
            )
        )

    if args.pad:
        logger.info('Padding positions')

        if args.interaction == 'peptide':
            logger.debug('Padding peptides')
            contact_map['class_i_pad_pos'] = None
            contact_map.loc[contact_map['relative_pos_centre_pmhc'].notna(), 'class_i_pad_pos'] = contact_map[
                contact_map['relative_pos_centre_pmhc'].notna()
            ]['relative_pos_centre_pmhc'].map(
                lambda pos: pos + args.peptide_length // 2 - 1,
            )

        else:
            logger.debug('Padding MHC')
            contact_map['imgt_pmhc'] = None
            shortened_resi = contact_map[contact_map['resi_pmhc'].notna()]['resi_pmhc'].map(
                lambda resi: (
                    (resi[0] - MHC_I_IMGT_BETA_HELIX_START, resi[1]) if resi[0] > MHC_I_IMGT_BETA_HELIX_START else resi
                )
            )
            mhc_helix = contact_map[contact_map['resi_pmhc'].notna()].apply(
                lambda row: assign_helix(
                    row.mhc_type, row.chain_type_pmhc, ''.join(str(element) for element in row.resi_pmhc)
                ),
                axis=1,
            )
            contact_map.loc[contact_map['resi_pmhc'].notna(), 'imgt_pmhc'] = pd.concat(
                [mhc_helix.to_frame(), shortened_resi.to_frame()],
                axis=1,
            ).apply(tuple, axis=1)

            logger.debug('Loading MHC pseudo sequence positions from %s', args.mhc_pseudo_seq_positions)
            with open(args.mhc_pseudo_seq_positions, 'r') as fh:
                mhc_pseudo_sequences = json.load(fh)

            mhc_pseudo_pos = [
                (
                    helix,
                    (
                        int(''.join([char for char in resi if char.isnumeric()])),
                        '' if resi[-1].isnumeric() else resi[-1],
                    ),
                )
                for helix in ('alpha', 'beta')
                for resi in mhc_pseudo_sequences[helix]
            ]

        logger.debug('Padding CDRs')
        contact_map['cdr_pad_pos'] = None
        for cdr_num, cdr_length in enumerate((args.cdr_1_length, args.cdr_2_length, args.cdr_3_length), 1):
            selection = (contact_map['cdr'] == cdr_num) & contact_map['relative_pos_centre_tcr'].notna()
            contact_map.loc[selection, 'cdr_pad_pos'] = contact_map[selection]['relative_pos_centre_tcr'].map(
                lambda pos, cdr_length=cdr_length: pos + cdr_length // 2 - 1,
            )

    overall_title = ['Peptides' if args.interaction == 'peptide' else 'MHC']

    overall_title.append(f"({'/'.join(args.mhc_types)})")

    if args.smooth:
        overall_title.append('- Gaussian')

    if args.pad:
        overall_title.append('- Padded')
        tcr_col_name = 'cdr_pad_pos'
        pmhc_col_name = 'class_i_pad_pos' if args.interaction == 'peptide' else 'imgt_pmhc'

    else:
        tcr_col_name = 'relative_pos_centre_tcr'
        pmhc_col_name = 'relative_pos_centre_pmhc' if args.interaction == 'peptide' else 'resi_pmhc'

    fig_data = contact_map[contact_map['mhc_type'].isin(args.mhc_types)]

    if args.interaction == 'peptide':
        fig_data = fig_data[fig_data['chain_type_pmhc'] == 'antigen_chain']

    else:
        fig_data = fig_data[fig_data['chain_type_pmhc'] != 'antigen_chain']

    fig_data = (
        fig_data.groupby(
            ['chain_type_tcr', 'cdr', tcr_col_name, pmhc_col_name],
        )['proportion']
        .agg('sum')
        .reset_index()
    )

    y_labels = sorted(fig_data[tcr_col_name].unique())
    x_labels = sorted(fig_data[pmhc_col_name].unique())

    tables = [[None, None, None], [None, None, None]]
    for i, chain_type in enumerate(('alpha', 'beta')):
        for j, (cdr, cdr_length) in enumerate(((1, args.cdr_1_length), (2, args.cdr_2_length), (3, args.cdr_3_length))):
            selected_data = fig_data[(fig_data['chain_type_tcr'] == chain_type) & (fig_data['cdr'] == cdr)]

            if args.pad:
                y_labels = np.arange(cdr_length)
                selected_data = selected_data[
                    (selected_data['cdr_pad_pos'] >= 0) & (selected_data['cdr_pad_pos'] < cdr_length)
                ]

                if args.interaction == 'peptide':
                    x_labels = np.arange(args.peptide_length)

                    selected_data = selected_data[
                        (selected_data['class_i_pad_pos'] >= 0)
                        & (selected_data['class_i_pad_pos'] < args.peptide_length)
                    ]

                else:
                    x_labels = mhc_pseudo_pos
                    selected_data = selected_data[selected_data['imgt_pmhc'].isin(mhc_pseudo_pos)]

            if args.smooth:
                group_max = selected_data['proportion'].max()

                x_dist = selected_data.groupby([pmhc_col_name])['proportion'].sum().reset_index()
                x_mu = np.average(x_dist[pmhc_col_name], weights=x_dist['proportion'])
                x_var = np.average((x_dist[pmhc_col_name] - x_mu) ** 2, weights=x_dist['proportion'])

                y_dist = selected_data.groupby([tcr_col_name])['proportion'].sum().reset_index()
                y_mu = np.average(y_dist[tcr_col_name], weights=y_dist['proportion'])
                y_var = np.average((y_dist[tcr_col_name] - y_mu) ** 2, weights=y_dist['proportion'])

                means = np.array([x_mu, y_mu])
                cov = np.array([[x_var, 0], [0, y_var]])

                gaussian = []
                for x in x_labels:
                    for y in y_labels:
                        gauss_xy = gaussian_2d(np.array([x, y]), means, cov)
                        gaussian.append([x, y, gauss_xy])

                selected_data = pd.DataFrame(gaussian, columns=[pmhc_col_name, tcr_col_name, 'proportion'])
                selected_data['proportion'] = selected_data['proportion'] * (
                    group_max / selected_data['proportion'].max()
                )

            table = selected_data.pivot_table(
                index=tcr_col_name,
                columns=pmhc_col_name,
                values='proportion',
            )

            table[[label for label in x_labels if label not in table.columns]] = np.nan
            missing_table = pd.DataFrame(
                np.nan,
                index=pd.Series([label for label in y_labels if label not in table.index], name=table.index.name),
                columns=table.columns,
            )
            if not missing_table.empty:
                table = pd.concat([table, missing_table])

            table = table[sorted(table.columns)].sort_index(ascending=False)
            table = table.fillna(0)

            if args.interaction == 'MHC':
                table.columns = (
                    [f'{seq_id}{insert_code} - {helix}' for helix, (seq_id, insert_code) in table.columns]
                    if args.pad
                    else [f'{seq_id}{insert_code}' for seq_id, insert_code in table.columns]
                )

            tables[i][j] = table

    if args.normalise:
        logger.debug('Normalising values')
        total = 0
        for i in range(2):
            for j in range(3):
                total += tables[i][j].sum().sum()

        for i in range(2):
            for j in range(3):
                tables[i][j] = tables[i][j].map(lambda val, total=total: val / total)

    if args.smooth or args.pad:
        maximum = 0
        for i in range(2):
            for j in range(3):
                maximum = max(maximum, tables[i][j].max().max())

    else:
        maximum = fig_data['proportion'].max()

    if args.separate_plots:
        for i, chain_type in enumerate(('alpha', 'beta')):
            for j, cdr in enumerate((1, 2, 3)):
                cdr_name = f'CDR{cdr}$\\{chain_type}$'

                plt.figure()
                ax = sns.heatmap(tables[i][j], vmax=maximum)
                ax.set_title(cdr_name)

                plt.tight_layout()

                if args.output:
                    dir_name = os.path.dirname(args.output)
                    base_name = os.path.basename(args.output)

                    plt.savefig(os.path.join(dir_name, cdr_name.replace('$', '').replace('\\', '') + '_' + base_name))
                    plt.close()

                else:
                    plt.show()

    else:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(' '.join(overall_title), fontsize=20)

        for i, chain_type in enumerate(('alpha', 'beta')):
            for j, cdr in enumerate((1, 2, 3)):
                sns.heatmap(tables[i][j], vmax=maximum, ax=axes[i, j])
                axes[i, j].set_title(f'CDR{cdr}$\\{chain_type}$')

        plt.tight_layout()

        if args.output:
            plt.savefig(args.output)

        else:
            plt.show()


if __name__ == '__main__':
    main()
