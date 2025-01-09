"""Create contact maps between TCRs and pMHCs."""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.chemistry import HEAVY_ATOMS
from tcr_antigen_prediction.data.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, assign_cdr_number
from tcr_antigen_prediction.data.structure import bio_to_pandas

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument('input', nargs='+', help='paths to the input PDB files')
inputs.add_argument('--summary-csv', required=True, help='path to the summary csv describing the PDB files')

output = parser.add_argument_group('Output')
output.add_argument('--output', '-o', required=True, help='path to output CSV')

parser.add_argument(
    '--cutoff-distance',
    default=5.0,
    type=float,
    help='maximum distance to consider a contact (default: 5.0 Å)',
)
parser.add_argument(
    '--tcr-norm',
    nargs='?',
    choices=[None, 'cdr_type', 'imgt_number', 'cdr_position', 'relative_pos_centre'],
    default='cdr_type',
    help="Normalisation strategy for the TCR counts (Default: 'cdr_type').",
)
parser.add_argument(
    '--mhc-norm',
    nargs='?',
    choices=[None, 'mhc_type', 'chain_type', 'imgt_number'],
    default='chain_type',
    help="Normalisation strategy for the MHC counts (Default: 'chain_type').",
)
parser.add_argument(
    '--peptide-norm',
    nargs='?',
    choices=[None, 'peptide_position', 'relative_pos_centre'],
    default='peptide_position',
    help="Normalisation strategy for the peptide counts (Default: 'peptide_position').",
)
parser.add_argument(
    '--even-offset-side',
    choices=['left', 'right'],
    default='left',
    help="Offset side when relative_pos_centre is selected (Default: 'left')",
)
parser.add_argument('--probability-filter', type=float, default=None, help='Filter out contacts below this cutoff')

add_logging_arguments(parser)


def enumerate_chain(chain_df: pd.DataFrame) -> pd.Series:
    """Enumerate a chain by residue to get a residue number."""
    chain_df = chain_df.copy()
    chain_df['residue_insert_code'] = chain_df['residue_insert_code'].fillna('')

    return chain_df.groupby(['residue_seq_id', 'residue_insert_code']).ngroup() + 1


def number_relative_to_centre(numbering: np.array, even_offset='left') -> np.array:
    """Numbers a sequential array relative to its centre.

    Examples:
        # Odd length example
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6, 7]))
        array([-3, -2, -1, 0, 1, 2, 3])

        # Even length example - left offset
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6]), even_offset='left')
        array([-2, -1, 0, 1, 2, 3])

        # Even length example - right offset
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6]), even_offset='left')
        array([-3, -2, -1, 0, 1, 2])

    Args:
        numbering: sequential input array
        even_offset: offset to give even length input numbering, can be either 'left' or 'right' (Default: 'left')

    Returns:
        array of the same length but with numbering relative to centre.

    Raises:
        ValueError: if even_offset is not 'left' or 'right'

    """
    length = np.max(numbering)

    if length % 2 == 0:
        match even_offset:
            case 'left':
                return numbering - length // 2

            case 'right':
                return numbering - (length // 2 + 1)

            case _:
                msg = f"Invalid even_offset={even_offset}. Must be 'left' or 'right'."
                raise ValueError(msg)

    return numbering - np.ceil(length / 2)


def collect_contacts(
    paths: list[str],
    summary_df: pd.DataFrame,
    cutoff_distance: float,
    even_offset_side: str = 'left',
) -> pd.DataFrame:
    """Collect the contacting residues from the structures."""
    contacts = []
    pdb_parser = PDBParser(QUIET=True)

    for path in paths:
        base_path = os.path.basename(path)
        logger.debug('Working on %s', base_path)

        meta_data = summary_df.loc[base_path]

        structure = pdb_parser.get_structure('', path)
        structure_df = bio_to_pandas(structure)

        logger.debug('Annotating structure')
        chain_annotations = {
            chain_id: chain_type
            for chain_type, chain_id in meta_data.filter(regex='.*chain[0-9]?').to_dict().items()
            if chain_id
        }

        structure_df['resi'] = structure_df['residue_seq_id'].apply(str) + structure_df['residue_insert_code'].fillna(
            ''
        )
        structure_df['chain_type'] = structure_df['chain_id'].map(chain_annotations)
        structure_df['cdr'] = None
        structure_df.loc[
            (structure_df['chain_type'] == 'Achain') | (structure_df['chain_type'] == 'Bchain'),
            'cdr',
        ] = structure_df[(structure_df['chain_type'] == 'Achain') | (structure_df['chain_type'] == 'Bchain')][
            'residue_seq_id'
        ].map(assign_cdr_number)

        logger.debug('Splitting TCR and pMHC')
        tcr_cdrs_df = structure_df[structure_df['cdr'].notna() & structure_df['element'].isin(HEAVY_ATOMS)].copy()

        tcr_cdrs_df['cdr_position'] = (
            tcr_cdrs_df.groupby(['chain_type', 'cdr'])
            .apply(
                enumerate_chain,
                include_groups=False,
            )
            .droplevel(['chain_type', 'cdr'])
        )

        tcr_cdrs_df['relative_pos_centre'] = tcr_cdrs_df.groupby(
            ['chain_type', 'cdr'],
        )['cdr_position'].transform(number_relative_to_centre, even_offset=even_offset_side)

        match meta_data['mhc_type']:
            case 'MH1':
                mhc_abd_df = structure_df[
                    (structure_df['chain_type'] == 'mhc_chain1')
                    & structure_df['residue_seq_id'].isin(IMGT_MH1_ABD)
                    & structure_df['element'].isin(HEAVY_ATOMS)
                ].copy()

            case 'MH2':
                mhc_abd_df = structure_df[
                    ((structure_df['chain_type'] == 'mhc_chain1') | (structure_df['chain_type'] == 'mhc_chain2'))
                    & structure_df['residue_seq_id'].isin(IMGT_MH2_ABD)
                    & structure_df['element'].isin(HEAVY_ATOMS)
                ].copy()

            case _:
                msg = f"Invalid MHC Type: {meta_data['mhc_type']}"
                raise ValueError(msg)

        peptide_df = structure_df[
            (structure_df['chain_type'] == 'antigen_chain') & structure_df['element'].isin(HEAVY_ATOMS)
        ].copy()

        peptide_df['peptide_position'] = enumerate_chain(peptide_df)
        peptide_df['relative_pos_centre'] = number_relative_to_centre(
            peptide_df['peptide_position'],
            even_offset=even_offset_side,
        )
        peptide_df['resi'] = None

        pmhc_df = pd.concat([peptide_df, mhc_abd_df])

        logger.debug('Calculating distances')
        interface = tcr_cdrs_df.merge(pmhc_df, how='cross', suffixes=('_tcr', '_pmhc'))
        interface['distance'] = np.sqrt(
            np.square(interface['pos_x_tcr'] - interface['pos_x_pmhc'])
            + np.square(interface['pos_y_tcr'] - interface['pos_y_pmhc'])
            + np.square(interface['pos_z_tcr'] - interface['pos_z_pmhc'])
        )

        logger.debug('Getting atoms within %.2f Å of each other', cutoff_distance)
        interface_contacts = interface[interface['distance'] < cutoff_distance].copy()
        interface_contacts['mhc_type'] = meta_data['mhc_type']
        interface_contacts['path'] = base_path
        contacts.append(interface_contacts)

    return pd.concat(contacts)


def get_normalisation_columns(
    tcr_strategy: str | None,
    mhc_strategy: str | None,
    peptide_strategy: str | None,
) -> list[str]:
    """Get normalisation columns based on the selected strategies."""
    norm_columns = []
    pmhc_selection_columns = []

    match tcr_strategy:
        case 'cdr_type':
            norm_columns.append('cdr_type')

        case 'imgt_number':
            norm_columns.append('cdr_type')
            norm_columns.append('resi_tcr')

        case 'cdr_position':
            norm_columns.append('cdr_type')
            norm_columns.append('cdr_position')

        case 'relative_pos_centre':
            norm_columns.append('cdr_type')
            norm_columns.append('relative_pos_centre_tcr')

    match mhc_strategy:
        case 'mhc_type':
            pmhc_selection_columns += ['mhc_chain1', 'mhc_chain2']
            norm_columns.append('mhc_type')

        case 'chain_type':
            pmhc_selection_columns += ['mhc_chain1', 'mhc_chain2']
            norm_columns.append('mhc_type')
            norm_columns.append('chain_type_pmhc')

        case 'imgt_number':
            pmhc_selection_columns += ['mhc_chain1', 'mhc_chain2']
            norm_columns.append('mhc_type')
            norm_columns.append('chain_type_pmhc')
            norm_columns.append('resi_pmhc')

    match peptide_strategy:
        case 'peptide_position':
            pmhc_selection_columns.append('antigen_chain')
            norm_columns.append('peptide_position')

        case 'relative_pos_centre':
            pmhc_selection_columns.append('antigen_chain')
            norm_columns.append('relative_pos_centre_pmhc')

    return norm_columns, pmhc_selection_columns


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    summary_df = pd.read_csv(args.summary_csv)

    base_paths = [os.path.basename(path) for path in args.input]
    summary_df = summary_df[summary_df['path'].isin(base_paths)].copy().set_index('path', drop=True)

    logger.info('Collecting contacts')
    contacts = collect_contacts(args.input, summary_df, args.cutoff_distance, even_offset_side=args.even_offset_side)

    logger.info('Normalising outputs')
    contacts = contacts.drop_duplicates(
        [
            'path',
            'chain_type_tcr',
            'residue_seq_id_tcr',
            'residue_insert_code_tcr',
            'chain_type_pmhc',
            'residue_seq_id_pmhc',
            'residue_insert_code_pmhc',
        ]
    )
    contacts['cdr_type'] = (
        'CDR'
        + contacts['cdr_tcr'].apply(int).apply(str)
        + contacts['chain_type_tcr'].map({'Achain': 'alpha', 'Bchain': 'beta'})
    )

    norm_columns, pmhc_selection_columns = get_normalisation_columns(args.tcr_norm, args.mhc_norm, args.peptide_norm)

    contacts = contacts[contacts['chain_type_pmhc'].isin(pmhc_selection_columns)]

    contact_counts = contacts.value_counts(norm_columns, dropna=False, normalize=True)

    if args.probability_filter:
        logger.info('Filtering probabilities under %.2f', args.probability_filter)
        contact_counts = contact_counts[contact_counts >= args.probability_filter]

    logger.info('Outputting results to %s', args.output)
    contact_counts = contact_counts.to_frame().reset_index()
    contact_counts = contact_counts.sort_values(norm_columns).reset_index(drop=True)

    contact_counts.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
