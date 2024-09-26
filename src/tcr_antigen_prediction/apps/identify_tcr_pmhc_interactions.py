'''Identify the interacting pairs of TCR:pMHC molecules in a PDB file.'''
import argparse
import logging
import re

import numpy as np
from Bio.PDB import PDBParser

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.imgt_numbering import IMGT_CDR, IMGT_VARIABLE_DOMAIN, IMGT_MH1_ABD, IMGT_MH2_ABD
from tcr_antigen_prediction.structure import get_header

logger = logging.getLogger()

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('structure', help='path to pdb structure')
parser.add_argument('--contact-distance', type=float, default=5.0,
                    help='distance to consider two chains part of the same complex (Default: 5.0)')
parser.add_argument('--output', '-o', help='path to output csv file')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    # Parse information
    with open(args.structure, 'r') as fh:
        header = get_header(fh.read())
        fh.seek(0)

        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('', fh)

    chain_maps = {chain_id: chain_type
                  for chain_type, chain_id in re.findall(r'(\w+)CHAINS?=(\w+)', header, flags=re.MULTILINE)}

    # Get relevant atom coordinates
    heavy_atom_coordinates = []
    chains = []
    residue_ids = []

    for model in structure:
        for chain in model:
            for res in chain:
                # Only consider TCR variable domains or MHC antigen binding domains
                if chain.id in chain_maps:
                    if ((chain_maps[chain.id] == 'A' or chain_maps[chain.id] == 'B')
                            and res.id[1] not in IMGT_VARIABLE_DOMAIN):
                        continue

                    if chain_maps[chain.id] == 'MH1' and res.id[1] not in IMGT_MH1_ABD:
                        continue

                    if (chain_maps[chain.id] == 'GA' or chain_maps[chain.id] == 'GB') and res.id[1] not in IMGT_MH2_ABD:
                        continue

                for atom in res:
                    if atom.element != 'H' and res.get_resname() != 'HOH':
                        heavy_atom_coordinates.append(atom.coord)
                        chains.append(chain.id)
                        residue_ids.append(res.id[1])

    heavy_atom_coordinates = np.vstack(heavy_atom_coordinates)
    chains = np.hstack(chains)
    residue_ids = np.hstack(residue_ids)

    # Annotate unknown chains
    chain_ids = np.unique(chains)

    for chain_id in chain_ids:
        if chain_id not in chain_maps:
            chain_maps[chain_id] = None

    # Find interacting chains
    interacting_chains = {}

    for chain_id in chain_ids:
        chain_idxs, = np.where(chains == chain_id)
        chain_coords = heavy_atom_coordinates[chain_idxs]

        other_chains_mask = np.ones(len(heavy_atom_coordinates), dtype=bool)
        other_chains_mask[chain_idxs] = False

        other_chain_coords = heavy_atom_coordinates[other_chains_mask]
        other_chain_ids = chains[other_chains_mask]

        distances = (
            np.sqrt(np.sum((chain_coords[:, np.newaxis, :] - other_chain_coords[np.newaxis, :, :]) ** 2, axis=2))
        )
        contacts = distances <= args.contact_distance
        contacting_other_coords = np.max(contacts, axis=0)

        contacting_other_chains = np.unique(other_chain_ids[contacting_other_coords]).tolist()

        interacting_chains[(chain_id,
                            chain_maps[chain_id])] = [(id_, chain_maps[id_]) for id_ in contacting_other_chains]

    # Separate complexes
    alpha_chains = [(chain_id, chain_type) for chain_id, chain_type in interacting_chains if chain_type == 'A']
    tcrs = []

    for chain_id, chain_type in alpha_chains:
        alpha_chain_interactions = interacting_chains[(chain_id, chain_type)]
        possible_beta_chains = [(chain_id, chain_type)
                                for chain_id, chain_type in alpha_chain_interactions
                                if chain_type == 'B']

        if len(possible_beta_chains) == 1:
            tcrs.append(((chain_id, chain_type), possible_beta_chains[0]))

        else:
            for option in possible_beta_chains:
                if (chain_id, chain_type) in interacting_chains[option]:
                    tcrs.append(((chain_id, chain_type), option))
                    break

            else:
                tcrs.append(((chain_id, chain_type), None))

    mhcs = []
    mhc_chain_1s = [(chain_id, chain_type) for chain_id, chain_type in interacting_chains if chain_type == 'MH1']

    if len(mhc_chain_1s) == 0:
        mhc_chain_1s = [(chain_id, chain_type) for chain_id, chain_type in interacting_chains if chain_type == 'GA']

    for chain_id, chain_type in mhc_chain_1s:
        mhc_chain2_type = 'B2M' if chain_type == 'MH1' else 'GB'
        chain_interactions = interacting_chains[(chain_id, chain_type)]
        possible_mhc_chain_2_interactions = [(chain_id, chain_type)
                                             for chain_id, chain_type in chain_interactions
                                             if chain_type == mhc_chain2_type]

        if len(possible_mhc_chain_2_interactions) == 1:
            mhcs.append(((chain_id, chain_type), possible_mhc_chain_2_interactions[0]))

        else:
            for option in possible_mhc_chain_2_interactions:
                if (chain_id, chain_type) in interacting_chains[option]:
                    mhcs.append(((chain_id, chain_type), option))
                    break

            else:
                mhcs.append(((chain_id, chain_type), None))

    # Add antigens
    unidentified_chains = [chain_id for chain_id, chain_type in chain_maps.items() if chain_type is None]
    pmhcs = []
    for mhc in mhcs:
        ((mhc_chain_1_id, mhc_chain_1_type), (mhc_chain_2_id, _)) = mhc

        relevant_mhc_chains = [mhc_chain_1_id] if mhc_chain_1_type == 'MH1' else [mhc_chain_1_id, mhc_chain_2_id]
        relevant_mhc_chains = np.array(relevant_mhc_chains)

        mhc_imgt_abd = IMGT_MH1_ABD if mhc_chain_1_type == 'MH1' else IMGT_MH2_ABD
        mhc_imgt_abd = np.array(sorted(list(mhc_imgt_abd)))

        mhc_mask = np.in1d(chains, relevant_mhc_chains)
        mhc_abd_mask = np.in1d(residue_ids, mhc_imgt_abd) & mhc_mask
        mhc_abd_coords = heavy_atom_coordinates[mhc_abd_mask]

        for chain_id in unidentified_chains:
            chain_mask = chains == chain_id
            chain_coords = heavy_atom_coordinates[chain_mask]

            distances = distances = (
                np.sqrt(np.sum((mhc_abd_coords[:, np.newaxis, :] - chain_coords[np.newaxis, :, :]) ** 2, axis=2))
            )
            contacts = distances <= args.contact_distance

            if contacts.max():
                pmhcs.append((mhc, (chain_id, 'AG')))
                break

    tcr_pmhcs = []
    for tcr, pmhc in ((tcr, pmhc) for tcr in tcrs for pmhc in pmhcs):
        (alpha_chain_id, _), (beta_chain_id, _) = tcr
        ((mhc_chain_1_id, mhc_chain_1_type), (mhc_chain_2_id, _)), (antigen_chain_id, _) = pmhc

        relevant_mhc_chains = [mhc_chain_1_id] if mhc_chain_1_type == 'MH1' else [mhc_chain_1_id, mhc_chain_2_id]
        relevant_mhc_chains = np.array(relevant_mhc_chains)

        mhc_imgt_abd = IMGT_MH1_ABD if mhc_chain_1_type == 'MH1' else IMGT_MH2_ABD
        mhc_imgt_abd = np.array(sorted(list(mhc_imgt_abd)))

        tcr_alpha_chain_mask = chains == alpha_chain_id
        tcr_beta_chain_mask = chains == beta_chain_id
        cdr_mask = (np.in1d(residue_ids, np.array(sorted(list(IMGT_CDR))))
                    & (tcr_alpha_chain_mask | tcr_beta_chain_mask))

        mhc_mask = np.in1d(chains, relevant_mhc_chains)
        mhc_abd_mask = np.in1d(residue_ids, mhc_imgt_abd) & mhc_mask

        antigen_mask = chains == antigen_chain_id

        cdr_coords = heavy_atom_coordinates[cdr_mask]
        antigen_coords = heavy_atom_coordinates[antigen_mask]

        distances = distances = (
            np.sqrt(np.sum((cdr_coords[:, np.newaxis, :] - antigen_coords[np.newaxis, :, :]) ** 2, axis=2))
        )

        contacts = distances <= args.contact_distance

        if contacts.max():
            tcr_pmhcs.append((tcr, pmhc))

    uncomplexed_tcrs = []
    for tcr in tcrs:
        for complexed_tcr, _ in tcr_pmhcs:
            if tcr == complexed_tcr:
                break
        else:
            uncomplexed_tcrs.append(tcr)

    uncomplexed_pmhcs = []
    for pmhc in pmhcs:
        for _, complexed_pmhc in tcr_pmhcs:
            if pmhc == complexed_pmhc:
                break
        else:
            uncomplexed_pmhcs.append(pmhc)

    output = [['Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2', 'mhc_type']]

    for (alpha_chain_id, _), (beta_chain_id, _) in uncomplexed_tcrs:
        output.append([alpha_chain_id, beta_chain_id, '', '', '', ''])

    for ((mhc_chain_1_id, _), (mhc_chain_2_id, _)), (antigen_chain_id, _) in uncomplexed_pmhcs:
        mhc_type = mhc_chain_1_id if mhc_chain_1_id == 'MH1' else 'MH2'
        output.append(['', '', antigen_chain_id, mhc_chain_1_id, mhc_chain_2_id, mhc_type])

    for tcr, pmhc in tcr_pmhcs:
        (alpha_chain_id, _), (beta_chain_id, _) = tcr
        ((mhc_chain_1_id, mhc_chain_1_type), (mhc_chain_2_id, _)), (antigen_chain_id, _) = pmhc
        mhc_type = mhc_chain_1_type if mhc_chain_1_type == 'MH1' else 'MH2'

        output.append([alpha_chain_id, beta_chain_id, antigen_chain_id, mhc_chain_1_id, mhc_chain_2_id, mhc_type])

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write('\n'.join([','.join(line) for line in output]))
            fh.write('\n')

    else:
        print('\n'.join(['\t'.join(line) for line in output]))


if __name__ == '__main__':
    main()
