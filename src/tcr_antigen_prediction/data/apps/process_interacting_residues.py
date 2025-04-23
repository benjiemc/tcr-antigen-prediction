"""Process interacting residue information into 20 x 20 matrice of amino acid contact probabilities."""

import argparse
import logging
import sys

import numpy as np
import pandas as pd
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.data.amino_acid_encodings import PROTEIN_LETTERS

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('contacts', help='path to input contact data file')
parser.add_argument(
    '--residue-column-names',
    nargs=2,
    required=True,
    help='Name of the two columns containing the resiudes formatted as three letter codes',
)
parser.add_argument(
    '--value-name',
    default='proportion',
    help="Name of the column containing the proportion made up of each residue pair (Default: 'proportion')",
)
parser.add_argument(
    '--selection-query',
    help=(
        'Custom pandas query to select certain values of the dataframe. See here for syntax: '
        'https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html.'
    ),
)
parser.add_argument('--output', '-o', required=True, help='File path to output matrix to (in txt form)')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    logger.info('Loading contact data')
    contacts = pd.read_csv(args.contacts)

    if args.selection_query:
        logger.info('Selecting values and re-normalising')
        contacts = contacts.query(args.selection_query).copy()
        contacts[args.value_name] = contacts[args.value_name] / contacts[args.value_name].sum()

    logger.info('Aggregating values')
    contacts = contacts.groupby(args.residue_column_names)[args.value_name].sum().reset_index()

    logger.info('Processing contact information')
    contacts[[col_name + '_olc' for col_name in args.residue_column_names]] = contacts[args.residue_column_names].map(
        lambda residue: IUPACData.protein_letters_3to1[residue.title()]
    )

    for col_name in args.residue_column_names:
        contacts[col_name + '_olc'] = pd.Categorical(
            contacts[col_name + '_olc'], categories=PROTEIN_LETTERS, ordered=True
        )

    contacts_table = contacts.pivot_table(
        index=args.residue_column_names[0] + '_olc',
        columns=args.residue_column_names[1] + '_olc',
        values=args.value_name,
        fill_value=0.0,
        observed=False,
    )

    logger.debug('Filling missing values')
    contacts_table = contacts_table.reindex(
        index=contacts[args.residue_column_names[0] + '_olc'].cat.categories,
        columns=contacts[args.residue_column_names[1] + '_olc'].cat.categories,
        fill_value=0.0,
    )

    logger.info('Outputting contacts...')
    np.savetxt(args.output, contacts_table.to_numpy())


if __name__ == '__main__':
    main()
