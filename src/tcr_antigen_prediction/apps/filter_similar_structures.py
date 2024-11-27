"""Filter TCR:pMHC structures with the same sequence (TCR CDRs and peptide) to an RMSD threshold."""

import argparse
import logging
import os
import sys

import pandas as pd
from Bio.PDB import PDBParser, Structure
from sklearn.cluster import AgglomerativeClustering

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.comparisons import compute_structural_distances

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structures', help='path to directory containing structures')
parser.add_argument('--output', '-o', help='path to output csv')
parser.add_argument('--summary-csv', required=True, help='path to summary csv file')
parser.add_argument(
    '--structural-similarity-cutoff',
    required=True,
    type=float,
    help='RMSD threshold for structures with the same CDR and peptide sequences (in Å)',
)

add_logging_arguments(parser)


def remove_similar_structures(
    df: pd.DataFrame, threshold: float, structures: dict[str, Structure.Structure]
) -> pd.DataFrame:
    """Remove structures with the same CDR and peptide sequences within the RMSD threshold."""
    output_dfs = []

    for (cdr_sequence, peptide_sequence, mhc_type), group in df.groupby(
        ['collated_cdrs', 'peptide_sequence', 'mhc_type']
    ):
        if len(group) == 1:
            output_dfs.append(group)
            continue

        logger.debug('Screening TCR: %s, peptide: %s, MHC: %s', cdr_sequence, peptide_sequence, mhc_type)

        sorted_group = group.sort_values('name')
        sorted_group = sorted_group.dropna(axis='columns', how='all')

        group_structures = [structures[name] for name in sorted_group['name'].to_numpy()]
        chain_maps = sorted_group.filter(regex=r'\w+chain\w*').to_dict('records')

        distance_matrix = compute_structural_distances(group_structures, chain_maps, mhc_type)
        clusters = (
            AgglomerativeClustering(
                metric='precomputed', distance_threshold=threshold, linkage='single', n_clusters=None
            )
            .fit(distance_matrix)
            .labels_
        )
        clusters = pd.Series(clusters, index=sorted_group.index)
        output_dfs.append(sorted_group[~clusters.duplicated()])

    return pd.concat(output_dfs)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    summary_df = pd.read_csv(args.summary_csv)
    summary_df['collated_cdrs'] = summary_df.filter(regex='cdr|CDR').apply('-'.join, axis=1)

    pdb_parser = PDBParser(QUIET=True)
    structures = {
        name: pdb_parser.get_structure(name, os.path.join(args.structures, name + '.pdb'))
        for name in summary_df['name'].tolist()
    }

    selected_structures = remove_similar_structures(summary_df, args.structural_similarity_cutoff, structures)

    output = open(args.output, 'w') if args.output else sys.stdout  # noqa: SIM115
    selected_structures.to_csv(output, index=False)

    if args.output:
        output.close()


if __name__ == '__main__':
    main()
