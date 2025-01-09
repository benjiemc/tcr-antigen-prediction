"""Calculate the TCR energy potentials for the TCRen model.

The original model was developped in this paper: https://www.nature.com/articles/s43588-024-00653-0.
"""

import argparse
import itertools
import logging
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import assign_cdr_number
from tcr_antigen_prediction.data.structure import bio_to_pandas

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('training_data', nargs='+', help='paths to training data for model')
parser.add_argument(
    '--strategy',
    choices=['regular', 'leave-one-out'],
    default='regular',
    help="training strategy employed to train the model (Default: 'regular')",
)
parser.add_argument(
    '--validation-data',
    nargs='+',
    help="path to validation data if validation is desired after fitting the model in 'regular' mode",
)
parser.add_argument('--summary-csv', required=True, help='path to summary csv file for the structures')
parser.add_argument('--output', '-o', required=True, help='path to output csv(s)')
parser.add_argument(
    '--contact-distance',
    default=5.0,
    type=float,
    help='threshold for contact distance between heavy atoms (Default: 5.0 Å)',
)

add_logging_arguments(parser)

AMINO_ACID_OLCS = list(IUPACData.protein_letters)


def load_data(paths: list[str], summary_df: pd.DataFrame, contact_distance: float) -> pd.DataFrame:
    """Get interacting residues from a list of TCR:pMHC PDB structures."""
    pdb_parser = PDBParser(QUIET=True)
    interacting_residues = defaultdict(list)

    for path in paths:
        entry_name = os.path.basename(path).replace('.pdb', '')
        logger.debug('Collecting contacts from %s', entry_name)

        pdb_id, chains = entry_name.split('_')

        row = summary_df[(summary_df['pdb'] == pdb_id) & (summary_df['chains'] == chains)].iloc[0]

        structure = pdb_parser.get_structure('', path)
        structure_df = bio_to_pandas(structure)

        chain_map = {val: key for key, val in row.filter(regex=r'\w+chain[1-2]?').to_dict().items()}
        structure_df['chain_type'] = structure_df['chain_id'].map(chain_map)
        structure_df['cdr'] = None
        structure_df['merge_key'] = 0

        tcr_selection = (structure_df['chain_type'] == 'Achain') | (structure_df['chain_type'] == 'Bchain')
        structure_df.loc[tcr_selection, 'cdr'] = structure_df[tcr_selection]['residue_seq_id'].map(assign_cdr_number)

        cdr_df = structure_df[tcr_selection & (structure_df['cdr'].notna())]
        antigen_df = structure_df[structure_df['chain_type'] == 'antigen_chain']

        # TODO replace with 'cross' when updating to a newer python version
        interaction = cdr_df.merge(antigen_df, on='merge_key', suffixes=('_cdr', '_antigen'))
        interaction = interaction.drop(columns=['merge_key'])

        coords_cdr = interaction.filter(regex=r'pos_[xyz]_cdr').to_numpy()
        coords_antigen = interaction.filter(regex=r'pos_[xyz]_antigen').to_numpy()

        interaction['distance'] = np.sqrt(np.sum((coords_cdr - coords_antigen) ** 2, axis=1))
        contacts = interaction[interaction['distance'] <= contact_distance]
        residue_contacts = contacts.drop_duplicates(
            [
                'chain_id_cdr',
                'residue_seq_id_cdr',
                'residue_insert_code_cdr',
                'chain_id_antigen',
                'residue_seq_id_antigen',
                'residue_insert_code_antigen',
            ]
        )

        interacting_residues['pdb'] += [row['pdb']] * len(residue_contacts)
        interacting_residues['Achain'] += [row['Achain']] * len(residue_contacts)
        interacting_residues['Bchain'] += [row['Bchain']] * len(residue_contacts)
        interacting_residues['antigen_chain'] += [row['antigen_chain']] * len(residue_contacts)
        interacting_residues['mhc_chain1'] += [row['mhc_chain1']] * len(residue_contacts)
        interacting_residues['mhc_chain2'] += [row['mhc_chain2']] * len(residue_contacts)

        interacting_residues['group_id'] += [row['group_id']] * len(residue_contacts)

        interacting_residues['tcr_cdr_residue'] += (
            residue_contacts['residue_name_cdr']
            .str.title()
            .map(
                IUPACData.protein_letters_3to1,
            )
            .tolist()
        )
        interacting_residues['peptide_residue'] += (
            residue_contacts['residue_name_antigen']
            .str.title()
            .map(
                IUPACData.protein_letters_3to1,
            )
            .tolist()
        )

    return pd.DataFrame.from_dict(interacting_residues)


def fit_tcr_en(interacting_residues: pd.DataFrame) -> pd.DataFrame:
    """Calculate the TCRen score from a table on interacting residues."""
    logger.debug('Calculating observed pairing probabilities')
    p_obs = interacting_residues.value_counts(normalize=True)
    p_obs.name = 'p_obs'

    logger.debug('Calculating expected pairing probabilities')
    p_a = interacting_residues.value_counts('tcr_cdr_residue', normalize=True)
    p_a.name = 'proportion'
    p_a = p_a.to_frame().reset_index()

    p_b = interacting_residues.value_counts('peptide_residue', normalize=True)
    p_b.name = 'proportion'
    p_b = p_b.to_frame().reset_index()

    # TODO replace with 'cross' when updating to a newer python version
    p_a['merge_key'] = 0
    p_b['merge_key'] = 0
    p_exp = p_a.merge(p_b, on='merge_key')
    p_exp = p_exp.drop(columns=['merge_key'])
    p_exp = p_exp.set_index(['tcr_cdr_residue', 'peptide_residue'])
    p_exp = p_exp['proportion_x'] * p_exp['proportion_y']
    p_exp.name = 'p_exp'

    potential_pairings = list(itertools.product(AMINO_ACID_OLCS, AMINO_ACID_OLCS))

    for pairing in potential_pairings:
        if pairing not in p_obs:
            logger.debug('Pairing %s not found in observed pairings. Assigning probability to 0.00', ':'.join(pairing))
            p_obs.loc[pairing] = 0.00

        if pairing not in p_exp:
            logger.debug('Pairing %s not found in expected pairings. Assigning probability to 0.00', ':'.join(pairing))
            p_exp.loc[pairing] = 0.00

    logger.debug('Calculating TCRen potentials')
    tcren = p_obs.to_frame().join(p_exp.to_frame())

    with np.errstate(divide='ignore'):
        tcren['tcren'] = -np.log(tcren['p_obs'] / tcren['p_exp'])

    tcren['tcren'] = tcren['tcren'].fillna(np.inf)

    return tcren[['tcren']]


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    summary_df = pd.read_csv(args.summary_csv)
    summary_df['chains'] = summary_df.apply(
        lambda row: ''.join([chain for chain in row.filter(regex=r'\w+chain[1-2]?').tolist() if not pd.isna(chain)]),
        axis=1,
    )

    if args.strategy == 'regular':
        logger.info('Loading data and finding contacting residues')
        training_data = load_data(args.training_data, summary_df, args.contact_distance)

        logger.info('Fitting potentials to training data')
        tcren = fit_tcr_en(training_data[['tcr_cdr_residue', 'peptide_residue']])

        if args.validation_data:
            logger.info('Evaluating TCRen model')
            logger.info('Loading data and finding contacting residues')
            evaluation_data = load_data(args.validation_data, summary_df, args.contact_distance)

            evaluation_tcr_ens = (
                evaluation_data.merge(
                    tcren.reset_index(),
                    how='left',
                )
                .groupby(
                    ['pdb', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2'],
                    dropna=False,
                )['tcren']
                .sum()
            )

            for (
                pdb_id,
                alpha_chain,
                beta_chain,
                antigen_chain,
                mhc_chain1,
                mhc_chain2,
            ), eval_tcren in evaluation_tcr_ens.items():
                logger.info(
                    (
                        'Evaluating '
                        'pdb=%s, Achain=%s, Bchain=%s, antigen_chain=%s, mhc_chain1=%s, mhc_chain2=%s: '
                        'TCRen=%f'
                    ),
                    pdb_id,
                    alpha_chain,
                    beta_chain,
                    antigen_chain,
                    mhc_chain1,
                    mhc_chain2,
                    eval_tcren,
                )

        logger.info('Outputting pottentials to %s', args.output)
        tcren.sort_values(['tcr_cdr_residue', 'peptide_residue']).to_csv(args.output)

    elif args.strategy == 'leave-one-out':
        logger.info('Loading data and finding contacting residues')
        summary_df['structure_name'] = summary_df['pdb'] + '_' + summary_df['chains']

        relevant_data_names = [os.path.basename(path).replace('.pdb', '') for path in args.training_data]
        relevant_summary_df = summary_df[summary_df['structure_name'].isin(relevant_data_names)]
        relevant_data = load_data(args.training_data, relevant_summary_df, args.contact_distance)

        groups = relevant_summary_df['group_id'].unique()

        for group in groups:
            logger.info('Fitting potentials to hold-out group: %d', group)
            training_data = relevant_data[relevant_data['group_id'] != group]
            evaluation_data = relevant_data[relevant_data['group_id'] == group]

            tcren = fit_tcr_en(training_data[['tcr_cdr_residue', 'peptide_residue']])

            logger.info('Evaluating leave-one-out group: %d', group)
            evaluation_tcr_ens = (
                evaluation_data.merge(
                    tcren.reset_index(),
                    how='left',
                )
                .groupby(
                    ['pdb', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2'],
                    dropna=False,
                )['tcren']
                .sum()
            )

            for (
                pdb_id,
                alpha_chain,
                beta_chain,
                antigen_chain,
                mhc_chain1,
                mhc_chain2,
            ), eval_tcren in evaluation_tcr_ens.items():
                logger.info(
                    (
                        'Evaluating group_id=%d, '
                        'pdb=%s, Achain=%s, Bchain=%s, antigen_chain=%s, mhc_chain1=%s, mhc_chain2=%s: '
                        'TCRen=%f'
                    ),
                    group,
                    pdb_id,
                    alpha_chain,
                    beta_chain,
                    antigen_chain,
                    mhc_chain1,
                    mhc_chain2,
                    eval_tcren,
                )

            base_output_name = os.path.basename(args.output)
            output_name = os.path.join(
                os.path.dirname(args.output),
                f"{base_output_name.split('.', 1)[0]}_LOO_{group}.{base_output_name.split('.', 1)[-1]}",
            )
            logger.info('Outputting potentials to %s', output_name)
            tcren.sort_values(['tcr_cdr_residue', 'peptide_residue']).to_csv(output_name)


if __name__ == '__main__':
    main()
