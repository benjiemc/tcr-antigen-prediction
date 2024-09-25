'''
Calculate the TCR energy potentials for the TCRen model. The original model was developped in this paper:
https://www.nature.com/articles/s43588-024-00653-0.
'''
import argparse
import itertools
import logging
import os
from typing import List

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.imgt_numbering import assign_cdr_number
from tcr_antigen_prediction.structure import bio_to_pandas

logger = logging.getLogger()

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('training_data', nargs='+', help='paths to training data for model')
parser.add_argument('--summary-csv', required=True, help='path to summary csv file for the structures')
parser.add_argument('--output', '-o', required=True, help='path to output csv')
parser.add_argument('--contact-distance', default=5.0, type=float,
                    help='threshold for contact distance between heavy atoms (Default: 5.0 Å)')

add_logging_arguments(parser)

AMINO_ACID_OLCS = list(IUPACData.protein_letters)


def load_data(paths: List[str], summary_df: pd.DataFrame, contact_distance: float) -> pd.DataFrame:
    '''Get interacting residues from a list of TCR:pMHC PDB structures.'''
    pdb_parser = PDBParser(QUIET=True)
    interacting_residues = []

    for path in paths:
        entry_name = os.path.basename(path).replace('.pdb', '')
        logger.info('Collecting contacts from %s', entry_name)

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

        cdr_df = structure_df[tcr_selection & (structure_df['cdr'].notnull())]
        antigen_df = structure_df[structure_df['chain_type'] == 'antigen_chain']

        # TODO replace with 'cross' when updating to a newer python version
        interaction = pd.merge(cdr_df, antigen_df, on='merge_key', suffixes=('_cdr', '_antigen'))
        interaction = interaction.drop(columns=['merge_key'])

        coords_cdr = interaction.filter(regex=r'pos_[xyz]_cdr').to_numpy()
        coords_antigen = interaction.filter(regex=r'pos_[xyz]_antigen').to_numpy()

        interaction['distance'] = np.sqrt(np.sum((coords_cdr - coords_antigen) ** 2, axis=1))
        contacts = interaction[interaction['distance'] <= contact_distance]
        residue_contacts = contacts.drop_duplicates([
            'chain_id_cdr', 'residue_seq_id_cdr', 'residue_insert_code_cdr',
            'chain_id_antigen', 'residue_seq_id_antigen', 'residue_insert_code_antigen',
        ])

        interacting_residues += list(zip(
            residue_contacts['residue_name_cdr'].str.title().map(IUPACData.protein_letters_3to1),
            residue_contacts['residue_name_antigen'].str.title().map(IUPACData.protein_letters_3to1),
        ))

    return pd.DataFrame(interacting_residues, columns=['tcr_cdr_residue', 'peptide_residue'])


def calculate_tcr_en(interacting_residues: pd.DataFrame) -> pd.DataFrame:
    '''Calculate the TCRen score from a table on interacting residues.'''
    logger.info('Calculating observed pairing probabilities')
    p_obs = interacting_residues.value_counts(normalize=True)
    p_obs.name = 'p_obs'

    logger.info('Calculating expected pairing probabilities')
    p_a = interacting_residues.value_counts('tcr_cdr_residue', normalize=True)
    p_a.name = 'proportion'
    p_a = p_a.to_frame().reset_index()

    p_b = interacting_residues.value_counts('peptide_residue', normalize=True)
    p_b.name = 'proportion'
    p_b = p_b.to_frame().reset_index()

    # TODO replace with 'cross' when updating to a newer python version
    p_a['merge_key'] = 0
    p_b['merge_key'] = 0
    p_exp = pd.merge(p_a, p_b, on='merge_key')
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

    logger.info('Calculating TCRen potentials')
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

    interacting_residues = load_data(args.training_data, summary_df, args.contact_distance)
    tcren = calculate_tcr_en(interacting_residues)

    logger.info('Outputting pottentials to %s', args.output)
    tcren.to_csv(args.output)


if __name__ == '__main__':
    main()
