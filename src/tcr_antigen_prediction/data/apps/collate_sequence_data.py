"""Collate TCR:pMHC sequence data from IEDB, McPAS-TCR, VDJdb, and ITRAP, removing redundant entries."""

import argparse
import json
import logging
import re
import sys
import warnings

import pandas as pd
import tidytcells
from Stitchr import stitchr as st
from Stitchr import stitchrfunctions as fxn

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3, MHC_I_IMGT_BETA_HELIX_START
from tcr_antigen_prediction.data.utils import mhc_slug_to_code

logger = logging.getLogger()

try:
    import anarci
except ImportError:
    logger.exception(
        'Some functions require ANARCI which is unavailable. see here for install instructions '
        'https://github.com/oxpig/ANARCI.'
    )

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('--output', '-o', help='path to the output csv file')

input_group = parser.add_argument_group('Input')
input_group.add_argument('--iedb-path', help='path to the IEDB data')
input_group.add_argument('--vdjdb-path', help='path to the VDJdb data')
input_group.add_argument('--itrap-path', help='path to the ITRAP data')
input_group.add_argument('--mcpas-tcr-path', help='path to the McPAS-TCR data')

annotations_group = parser.add_argument_group('Annotations')
annotations_group.add_argument(
    '--mhc-sequences',
    nargs='+',
    required=True,
    help='paths to mhc sequence information to create pseudo sequences.',
)
annotations_group.add_argument(
    '--mhc-pseudo-sequence-imgt-numbers',
    required=True,
    help='path to mhc pseudo sequence imgt numbers (JSON format)',
)

add_logging_arguments(parser)

COMMON_COLUMNS = [
    'cdr3_alpha',
    'v_alpha',
    'j_alpha',
    'cdr3_beta',
    'v_beta',
    'j_beta',
    'mhc_type',
    'mhc1',
    'mhc2',
    'peptide',
    'species',
]

REFERENCE_10X = (
    'https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-'
    'highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#'
)


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
    def create_input_args(args, gene_types):
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


def get_pseudo_sequence(
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


def process_iedb(iedb: pd.DataFrame) -> pd.DataFrame:
    """Process data from downloaded IEDB dataset."""
    logger.info('Number of sequences in IEDB: %d', len(iedb))

    logger.info('Selecting sequences with complete gene and CDR information')
    iedb = iedb[
        iedb['Chain 1 - CDR3 Calculated'].notna()
        & iedb['Chain 1 - Calculated V Gene'].notna()
        & iedb['Chain 1 - Calculated J Gene'].notna()
        & iedb['Chain 2 - CDR3 Calculated'].notna()
        & iedb['Chain 2 - Calculated V Gene'].notna()
        & iedb['Chain 2 - Calculated J Gene'].notna()
        & iedb['Epitope - Name'].notna()
        & iedb['Assay - MHC Allele Names'].notna()
        & (iedb['Receptor - Type'] == 'alphabeta')
    ]
    logger.info('Number of complete sequences: %d', len(iedb))

    logger.info('Selecting only peptide antigens')
    iedb = iedb[iedb['Epitope - Name'].str.contains(r'^[A-Z]+$', regex=True)]
    logger.info('Number of sequences with peptide antigens: %d', len(iedb))

    logger.info('Expanding MHC alleles')
    iedb['mhc_processed'] = iedb['Assay - MHC Allele Names'].str.split(', ')
    iedb = iedb.explode('mhc_processed')
    logger.info('Number of sequences after expanding MHC alleles: %d', len(iedb))

    logger.info('Splitting MH1 and MH2 gene names')
    iedb['mhc_processed'] = iedb['mhc_processed'].str.split(' ').map(lambda mhc: mhc[0])
    iedb[['mhc1', 'mhc2']] = iedb['mhc_processed'].str.split('/').apply(pd.Series)
    iedb['mhc_type'] = iedb['mhc1'].map(assign_mhc_class)
    iedb.loc[(iedb['mhc_type'] == 'MH1') & ~iedb['mhc2'].notna(), 'mhc2'] = 'B2M'
    # TODO check if this is ok:
    iedb.loc[(iedb['mhc_type'] == 'MH2') & ~iedb['mhc2'].notna(), 'mhc2'] = iedb.loc[
        (iedb['mhc_type'] == 'MH2') & ~iedb['mhc2'].notna()
    ]['mhc1']

    iedb = iedb.dropna(subset=['mhc_type', 'mhc1', 'mhc2'])
    logger.info('Number of sequences after splitting MH1 and MH2 gene names: %d', len(iedb))

    logger.info('Assigning species based on MHC gene name')
    iedb['species'] = iedb['mhc1'].map(assign_species)

    logger.info('Standardising column names')
    iedb = iedb.rename(
        {
            'Chain 1 - CDR3 Calculated': 'cdr3_alpha',
            'Chain 1 - Calculated V Gene': 'v_alpha',
            'Chain 1 - Calculated J Gene': 'j_alpha',
            'Chain 2 - CDR3 Calculated': 'cdr3_beta',
            'Chain 2 - Calculated V Gene': 'v_beta',
            'Chain 2 - Calculated J Gene': 'j_beta',
            'Epitope - Name': 'peptide',
        },
        axis='columns',
    )

    iedb = iedb[COMMON_COLUMNS].copy()
    iedb['source'] = 'IEDB'

    logger.info('Number of IEDB sequences: %d', len(iedb))

    return iedb


def process_vdjdb(vdjdb: pd.DataFrame) -> pd.DataFrame:
    """Process data from downloaded VDJdb dataset."""
    logger.info('Number of sequences in VDJdb: %d', len(vdjdb))

    logger.info('Merging paired alpha and beta TCR sequences')
    vdjdb_paired_alpha = vdjdb[(vdjdb['complex.id'] != 0) & (vdjdb['Gene'] == 'TRA')]
    vdjdb_paired_beta = vdjdb[(vdjdb['complex.id'] != 0) & (vdjdb['Gene'] == 'TRB')]

    logger.debug('Number of paired alpha-chain sequences: %d', len(vdjdb_paired_alpha))
    logger.debug('Number of paired beta-chain sequences: %d', len(vdjdb_paired_beta))

    vdjdb = vdjdb_paired_alpha.merge(
        vdjdb_paired_beta,
        how='inner',
        on='complex.id',
        suffixes=('_alpha', '_beta'),
    )

    duplicate_columns = [
        col
        for col in vdjdb_paired_alpha.columns.intersection(vdjdb_paired_beta.columns)
        if col != 'complex.id' and all(vdjdb[col + '_alpha'] == vdjdb[col + '_beta'])
    ]

    for col in duplicate_columns:
        vdjdb[col] = vdjdb[col + '_alpha']
        vdjdb = vdjdb.drop(columns=[col + '_alpha', col + '_beta'])

    logger.info('Number of paired alpha-beta TCR sequences: %d', len(vdjdb))

    logger.info('Removing 10x study')
    vdjdb = vdjdb[(vdjdb['Reference_alpha'] != REFERENCE_10X) & (vdjdb['Reference_beta'] != REFERENCE_10X)]

    logger.info('Number of sequences after removing 10x study: %d', len(vdjdb))

    logger.info('Standardising column names')
    vdjdb = vdjdb.rename(
        {
            'CDR3_alpha': 'cdr3_alpha',
            'V_alpha': 'v_alpha',
            'J_alpha': 'j_alpha',
            'CDR3_beta': 'cdr3_beta',
            'V_beta': 'v_beta',
            'J_beta': 'j_beta',
            'MHC class': 'mhc_type',
            'Epitope': 'peptide',
            'Species': 'species',
            'MHC A': 'mhc1',
            'MHC B': 'mhc2',
        },
        axis='columns',
    )

    vdjdb = vdjdb[COMMON_COLUMNS].copy()
    vdjdb['source'] = 'VDJdb'

    logger.info('Number of VDJdb sequences: %d', len(vdjdb))

    return vdjdb


def process_itrap(itrap: pd.DataFrame) -> pd.DataFrame:
    """Process data from downloaded ITRAP dataset."""
    logger.info('Number of sequences in ITRAP dataset: %d', len(itrap))

    logger.info('Standardising column names')
    itrap[['v_alpha', 'j_alpha', 'constant_alpha']] = itrap['genes_TRA'].str.split(';').apply(pd.Series)
    itrap[['v_beta', 'd_beta', 'j_beta', 'constant_beta']] = itrap['genes_TRB'].str.split(';').apply(pd.Series)

    itrap[['peptide', 'mhc1']] = itrap['peptide_HLA'].str.split().apply(pd.Series)

    itrap['mhc2'] = 'B2M'
    itrap['mhc_type'] = 'MH1'

    itrap['species'] = 'Human'

    itrap_renamed = itrap.rename(
        {
            'cdr3_TRA': 'cdr3_alpha',
            'cdr3_TRB': 'cdr3_beta',
        },
        axis='columns',
    )

    itrap_selected = itrap_renamed[COMMON_COLUMNS].copy()
    itrap_selected['source'] = 'ITRAP'

    logger.info('Number of ITRAP sequences: %d', len(itrap))

    return itrap_selected


def process_mcpas_tcr(mcpas_tcr: pd.DataFrame) -> pd.DataFrame:
    """Process data from downloaded McPAS-TCR dataset."""
    logger.info('Number of sequences in McPAS-TCR dataset: %d', len(mcpas_tcr))

    logger.info('Selecting sequences with complete gene and CDR information')
    mcpas_tcr = mcpas_tcr[
        mcpas_tcr['CDR3.alpha.aa'].notna()
        & mcpas_tcr['CDR3.beta.aa'].notna()
        & mcpas_tcr['Epitope.peptide'].notna()
        & mcpas_tcr['MHC'].notna()
        & mcpas_tcr['TRAV'].notna()
        & mcpas_tcr['TRAJ'].notna()
        & mcpas_tcr['TRBV'].notna()
        & mcpas_tcr['TRBJ'].notna()
    ].copy()

    logger.info('Number of sequences with complete gene and CDR information: %d', len(mcpas_tcr))

    logger.info('Standardising MHC gene names')
    mcpas_tcr['mhc_processed'] = mcpas_tcr['MHC'].str.replace(' ', '')
    mcpas_tcr['mhc_processed'] = mcpas_tcr['mhc_processed'].replace(
        {
            'HLA-A2:01': 'HLA-A*02:01',
            'HLA-A*2:01': 'HLA-A*02:01',
            'HLA-A2': 'HLA-A*02',
            'DRB1*04:01': 'HLA-DRB1*04:01',
            'HLA-A*011': 'HLA-A*11',
            'DQ8-trans': 'HLA-DQ8',
        }
    )

    mcpas_tcr['mhc_processed'] = mcpas_tcr['mhc_processed'].str.replace(r'^H-2', 'H2-', regex=True)
    mcpas_tcr['mhc2'] = 'B2M'

    logger.info('Assigning MHC type based on MHC allele codes')
    mcpas_tcr['mhc_type'] = mcpas_tcr['mhc_processed'].map(assign_mhc_class)

    logger.info('Standardising column names')
    mcpas_tcr = mcpas_tcr.rename(
        {
            'CDR3.alpha.aa': 'cdr3_alpha',
            'CDR3.beta.aa': 'cdr3_beta',
            'TRAV': 'v_alpha',
            'TRAJ': 'j_alpha',
            'TRBV': 'v_beta',
            'TRBJ': 'j_beta',
            'Epitope.peptide': 'peptide',
            'Species': 'species',
            'mhc_processed': 'mhc1',
        },
        axis='columns',
    )

    mcpas_tcr = mcpas_tcr[COMMON_COLUMNS].copy()
    mcpas_tcr['source'] = 'McPAS-TCR'

    logger.info('Number of McPAS-TCR Sequences: %d', len(mcpas_tcr))

    return mcpas_tcr


def collate_sequence_data(
    *datasets: pd.DataFrame,
    mhc_sequences: pd.DataFrame,
    mhc_pseudo_seq_imgt_positions: dict[str, list[int]],
) -> pd.DataFrame:
    """Collate datasets, standardise nomenclature, and assign sequences.

    Args:
        *datasets: datasets to collate together with the following columns
            - cdr3_alpha
            - v_alpha
            - j_alpha
            - cdr3_beta
            - v_beta
            - j_beta
            - mhc_type
            - mhc1
            - mhc2
            - peptide
            - species
        mhc_sequences: dataframe with MHC sequences indexed by MHC slug (simplified allele code)
        mhc_pseudo_seq_imgt_positions: dictionary with the imgt numbers for the alpha and beta helix pseudo sequence
            positions

    Returns:
        dataframe with the collated and standardised dataset

    """
    sequence_data = pd.concat(datasets).reset_index(drop=True)
    logger.info('Number of collated sequences: %d', len(sequence_data))

    logger.info('Standardising gene names')
    logger.debug('Standardising species names (Human/Mouse)')
    sequence_data['species'] = sequence_data['species'].replace(
        {
            'MusMusculus': 'Mouse',
            'HomoSapiens': 'Human',
        }
    )

    logger.debug('Standardising MHC types (MH1/MH2)')
    sequence_data['mhc_type'] = sequence_data['mhc_type'].replace(
        {'MHCI': 'MH1', 'MHCII': 'MH2'},
    )

    logger.debug('Standardising TCR gene names')
    sequence_data.loc[sequence_data['species'] == 'Human', ['v_alpha', 'j_alpha', 'v_beta', 'j_beta']] = (
        sequence_data.loc[sequence_data['species'] == 'Human', ['v_alpha', 'j_alpha', 'v_beta', 'j_beta']].map(
            tidytcells.tr.standardise,
            species='homosapiens',
            log_failures=False,
        )
    )

    sequence_data.loc[sequence_data['species'] == 'Mouse', ['v_alpha', 'j_alpha', 'v_beta', 'j_beta']] = (
        sequence_data.loc[sequence_data['species'] == 'Mouse', ['v_alpha', 'j_alpha', 'v_beta', 'j_beta']].map(
            tidytcells.tr.standardise,
            species='musmusculus',
            log_failures=False,
        )
    )

    sequence_data[['v_alpha', 'j_alpha', 'v_beta', 'j_beta']] = sequence_data[
        ['v_alpha', 'j_alpha', 'v_beta', 'j_beta']
    ].map(lambda gene: gene.split('/')[0] if gene else gene)

    logger.debug('Standardising CDR3 junctions')
    sequence_data[['cdr3_alpha', 'cdr3_beta']] = sequence_data[['cdr3_alpha', 'cdr3_beta']].map(
        tidytcells.junction.standardise,
        log_failures=False,
    )

    logger.debug('Standardising MHC gene names')
    sequence_data.loc[sequence_data['species'] == 'Human', ['mhc1', 'mhc2']] = sequence_data.loc[
        sequence_data['species'] == 'Human', ['mhc1', 'mhc2']
    ].map(
        tidytcells.mh.standardise,
        species='homosapiens',
        precision='protein',
        log_failures=False,
    )

    sequence_data.loc[sequence_data['species'] == 'Mouse', ['mhc1', 'mhc2']] = sequence_data.loc[
        sequence_data['species'] == 'Mouse', ['mhc1', 'mhc2']
    ].map(
        tidytcells.mh.standardise,
        species='musmusculus',
        precision='protein',
        log_failures=False,
    )

    sequence_data.loc[(sequence_data['mhc_type'] == 'MH1') & ~sequence_data['mhc2'].notna(), 'mhc2'] = 'B2M'

    sequence_data = sequence_data.dropna()
    logger.info('Number of sequences after standardising gene names (removing invalid names): %d', len(sequence_data))

    logger.info('Creating full length sequences from genes and junctions')
    sequence_data['alpha_chain_sequence'] = sequence_data.apply(
        lambda row: stitch_sequence(row.v_alpha, row.j_alpha, row.cdr3_alpha, row.species), axis=1
    )

    sequence_data['beta_chain_sequence'] = sequence_data.apply(
        lambda row: stitch_sequence(row.v_beta, row.j_beta, row.cdr3_beta, row.species), axis=1
    )

    sequence_data = sequence_data.dropna()
    logger.info('Number of sequences after creating full length sequences: %d', len(sequence_data))

    logger.info('IMGT numbering TCR sequences and extracting CDR regions')
    logger.debug('Numbering alpha-chains')
    sequence_data[['cdr1_alpha', 'cdr2_alpha', 'cdr3_alpha']] = (
        sequence_data['alpha_chain_sequence'].map(get_cdr_sequences).apply(pd.Series)
    )

    logger.debug('Numbering beta-chains')
    sequence_data[['cdr1_beta', 'cdr2_beta', 'cdr3_beta']] = (
        sequence_data['beta_chain_sequence'].map(get_cdr_sequences).apply(pd.Series)
    )

    sequence_data = sequence_data.dropna()
    logger.info('Number of sequences after extracting CDR regions: %d', len(sequence_data))

    logger.info('IMGT numbering MHC sequencs and extracting pseudo sequence representations')
    logger.debug('Adding MHC sequences')
    sequence_data = (
        sequence_data.merge(
            mhc_sequences[['mhc_code', 'canonical_sequence']],
            how='left',
            left_on='mhc1',
            right_on='mhc_code',
        )
        .rename({'canonical_sequence': 'mhc1_sequence'}, axis='columns')
        .merge(
            mhc_sequences[['mhc_code', 'canonical_sequence']],
            how='left',
            left_on='mhc2',
            right_on='mhc_code',
        )
        .rename({'canonical_sequence': 'mhc2_sequence'}, axis='columns')
    )

    sequence_data = sequence_data[sequence_data['mhc1_sequence'].notna()]
    sequence_data = sequence_data[sequence_data['mhc2_sequence'].notna() | (sequence_data['mhc_type'] == 'MH1')]

    logger.debug('Shortening to pseudo sequence')
    sequence_data['mhc_pseudo'] = sequence_data.apply(
        lambda row, mhc_pseudo_seq_imgt_positions=mhc_pseudo_seq_imgt_positions: get_pseudo_sequence(
            row.mhc1_sequence,
            row.mhc2_sequence,
            row.mhc_type,
            mhc_pseudo_seq_imgt_positions,
        ),
        axis='columns',
    )

    logger.info('Number of sequences after extracting MHC pseudo sequences: %d', len(sequence_data))

    logger.info('Collapsing redundant sequences')
    sequence_data = (
        sequence_data.groupby(
            [
                'cdr1_alpha',
                'cdr2_alpha',
                'cdr3_alpha',
                'cdr1_beta',
                'cdr2_beta',
                'cdr3_beta',
                'peptide',
                'mhc_pseudo',
            ],
            dropna=False,
        )[['v_alpha', 'j_alpha', 'v_beta', 'j_beta', 'mhc1', 'mhc2', 'mhc_type', 'species', 'source']]
        .agg(lambda sources: ';'.join(sorted(set(sources))))
        .reset_index()
    )
    logger.info('Number of non-redundant sequences: %d', len(sequence_data))

    return sequence_data


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.debug('Loading MHC Pseudo sequence IMGT positions')
    with open(args.mhc_pseudo_sequence_imgt_numbers, 'r') as fh:
        mhc_pseudo_seq_imgt_positions = json.load(fh)

    logger.debug('Loading MHC canonical sequences')
    mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in args.mhc_sequences])

    mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
    mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

    mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Human', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='homosapiens')

    mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Mouse', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='musmusculus')

    logger.info('Processing IEDB data')
    iedb = process_iedb(pd.read_csv(args.iedb_path))

    logger.info('Processing VDJdb data')
    vdjdb = process_vdjdb(pd.read_csv(args.vdjdb_path, delimiter='\t'))

    logger.info('Processing ITRAP data')
    itrap = process_itrap(pd.read_csv(args.itrap_path))

    logger.info('Processing McPAS-TCR data')
    mcpas_tcr = process_mcpas_tcr(pd.read_csv(args.mcpas_tcr_path))

    logger.info('Collating datasets')
    dataset = collate_sequence_data(
        iedb,
        vdjdb,
        itrap,
        mcpas_tcr,
        mhc_sequences=mhc_sequences,
        mhc_pseudo_seq_imgt_positions=mhc_pseudo_seq_imgt_positions,
    )

    logger.info('Outputting data to %s', args.output)
    dataset.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
