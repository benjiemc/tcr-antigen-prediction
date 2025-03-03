"""Annotate summary file with MHC Pseudo sequence information."""

import argparse
import json
import logging
import sys

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.imgt_numbering import MHC_I_IMGT_BETA_HELIX_START

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structures', nargs='+', help='')
parser.add_argument('--summary-csv', required=True, help='')
parser.add_argument('--mhc-pseudo-sequence-imgt-numbers', required=True, help='')
parser.add_argument('--output', '-o', required=True, help='')

add_logging_arguments(parser)


def add_mhc_pseudo_sequences(
    full_path: str,
    mhc_type: str,
    mhc_chain1: str | None,
    mhc_chain2: str | None,
    mhc_pseudo_imgt_numbers: dict[str, list[tuple[int, str]]],
) -> str:
    """Add MHC-TCR contact pseudo sequences."""
    logger.debug('Working on %s', full_path)

    structure = PDBParser().get_structure('', full_path)
    numbering = {
        chain_id: [
            (
                (res.id[1], res.id[2].strip()),
                IUPACData.protein_letters_3to1[res.get_resname().title()],
            )
            for res in structure[0][chain_id]
            if res.id[0] == ' '
        ]
        for chain_id in (mhc_chain1, mhc_chain2)
        if not pd.isna(chain_id)
    }

    mhc_pseudo_seq = []
    for helix, residues in mhc_pseudo_imgt_numbers.items():
        for imgt_seq_id, imgt_insert_code in residues:
            if helix == 'alpha':
                for (seq_id, insert_code), res_olc in numbering[mhc_chain1]:
                    if seq_id == imgt_seq_id and insert_code.strip() == imgt_insert_code:
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

            elif helix == 'beta' and mhc_type == 'MH1':
                for (seq_id, insert_code), res_olc in numbering[mhc_chain1]:
                    if (
                        seq_id == (imgt_seq_id + MHC_I_IMGT_BETA_HELIX_START)
                        and insert_code.strip() == imgt_insert_code
                    ):
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

            elif helix == 'beta' and mhc_type == 'MH2':
                for (seq_id, insert_code), res_olc in numbering[mhc_chain2]:
                    if seq_id == imgt_seq_id and insert_code.strip() == imgt_insert_code:
                        mhc_pseudo_seq.append(res_olc)
                        break

                else:
                    mhc_pseudo_seq.append('-')

    return ''.join(mhc_pseudo_seq)


def complete_path(name: str, paths: list[str]) -> str | None:
    full_paths = [path for path in paths if name in path]

    if len(full_paths) == 0:
        msg = f'No path given for {name}, defaulting to None'
        logger.warning(msg)
        return None

    return full_paths[0]


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading summary file')
    summary_df = pd.read_csv(args.summary_csv)

    logger.info('Finding paths for structures')
    summary_df['full_path'] = summary_df['path'].apply(complete_path, paths=args.structures)

    logger.info('Loading MHC Pseudo IMGT numbers')
    with open(args.mhc_pseudo_sequence_imgt_numbers, 'r') as fh:
        mhc_pseudo_imgt_numbers = {
            helix: [
                (
                    int(''.join([char for char in resi if char.isnumeric()])),
                    ''.join([char for char in resi if not char.isnumeric()]),
                )
                for resi in residues
            ]
            for helix, residues in json.load(fh).items()
        }

    logger.info('Collecting MHC Pseudo sequences')

    summary_df['mhc_pseudo'] = None
    summary_df.loc[summary_df['full_path'].notna(), 'mhc_pseudo'] = summary_df.loc[
        summary_df['full_path'].notna()
    ].apply(
        lambda row: add_mhc_pseudo_sequences(
            row.full_path,
            row.mhc_type,
            row.mhc_chain1,
            row.mhc_chain2,
            mhc_pseudo_imgt_numbers,
        ),
        axis=1,
    )

    logger.info('Outputting summary csv to %s', args.output)
    summary_df[[column for column in summary_df.columns if column != 'full_path']].to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
