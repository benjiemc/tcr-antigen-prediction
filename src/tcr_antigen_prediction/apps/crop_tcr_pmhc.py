"""Crop TCR-pMHC to binding interface (TCR Variable domain and MHC antigen binding domain)."""

import argparse
import logging
import sys

from Bio.PDB import PDBIO, Chain, Model, PDBParser, Select, Structure

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, IMGT_VARIABLE_DOMAIN

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to the pdb')
parser.add_argument('--output', '-o', required=True, help='path to the output pdb')
parser.add_argument('--tcr-chains', nargs='*', help='pdb chains for the tcr structure')
parser.add_argument('--mhc-chains', nargs='*', help='pdb chains for the mhc structure')
parser.add_argument('--antigen-chain', nargs='?', help='pdb chain for the antigen')
parser.add_argument('--remove-het-atoms', action='store_true', help='remove hetero atoms from structure')

add_logging_arguments(parser)


def crop_chain(chain: Chain.Chain, numbering: set[int]) -> Chain.Chain:
    """Crop chain based on numbering."""
    new_chain = Chain.Chain(chain.id)

    for residue in chain:
        if residue.id[1] in numbering or residue.id[0] != ' ':
            new_chain.add(residue.copy())

    return new_chain


class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == ' ' else 0


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    tcr_chains = args.tcr_chains if args.tcr_chains else []
    mhc_chains = args.mhc_chains if args.mhc_chains else []
    antigen_chain = args.antigen_chain if args.antigen_chain else None

    logger.info('Loading structure')
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('', args.structure)

    crop_structure = Structure.Structure('cropped')

    for model in structure:
        new_model = Model.Model(model.id)

        for chain in model:
            if chain.id in tcr_chains:
                logger.info('Cropping TCR chain %s', chain.id)
                new_chain = crop_chain(chain, IMGT_VARIABLE_DOMAIN)

            elif chain.id in mhc_chains:
                logger.info('Cropping MHC chain %s', chain.id)
                new_chain = crop_chain(chain, IMGT_MH1_ABD if len(args.mhc_chains) == 1 else IMGT_MH2_ABD)

            elif chain.id == antigen_chain and antigen_chain is not None:
                new_chain = chain.copy()

            else:
                logger.info('Removing chain %s', chain.id)
                continue

            new_model.add(new_chain)

        crop_structure.add(new_model)

    logger.info('Saving structure')
    io = PDBIO()
    io.set_structure(crop_structure)

    if args.remove_het_atoms:
        logger.info('Removing hetero atoms')
        logger.info('Saving structure')
        io.save(args.output, NonHetSelect())

    else:
        logger.info('Saving structure')
        io.save(args.output)


if __name__ == '__main__':
    main()
