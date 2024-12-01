"""Constants and functions for annotating sequences as CDR domains in T cell receptors."""

import logging

from Bio.PDB import Chain
from Bio.SeqUtils import IUPACData

logger = logging.getLogger(__name__)

try:
    import anarci
except ImportError:
    logger.exception(
        'Some functions require ANARCI which is unavailable. see here for install instructions '
        'https://github.com/oxpig/ANARCI).'
    )


IMGT_CDR1: set[int] = set(range(27, 38 + 1))
'''IMGT residue numbers corresponding to CDR 1 domains.'''
IMGT_CDR2: set[int] = set(range(56, 65 + 1))
'''IMGT residue numbers corresponding to CDR 2 domains.'''
IMGT_CDR3: set[int] = set(range(105, 117 + 1))
'''IMGT residue numbers corresponding to CDR 3 domains.'''
IMGT_CDR: set[int] = IMGT_CDR1.union(IMGT_CDR2).union(IMGT_CDR3)
'''IMGT residue numbers corresponding to all CDR domains.'''

IMGT_VARIABLE_DOMAIN: set[int] = set(range(1, 128 + 1))
'''Variable domain range for IMGT numbered TCR structures.'''

IMGT_FRAMEWORK_REGION = IMGT_VARIABLE_DOMAIN - IMGT_CDR
'''Framework (Fw) region range for IMGT numbered TCR structures.'''

IMGT_MH1_ABD: set[int] = set(range(1, 92)) | set(range(1001, 1092))
'''IMGT ranges of the antigen binding domain of MHC class I molecules.'''

IMGT_MH2_ABD: set[int] = set(range(1, 92))
'''IMGT ranges of the antigen binding domain of MHC class II molecules.'''


def assign_cdr_number(seq_id: int) -> int | None:
    """Assign CDR number for a sequence ID or return None if the sequence ID if the ID is not in CDR range."""
    if seq_id in IMGT_CDR1:
        return 1

    if seq_id in IMGT_CDR2:
        return 2

    if seq_id in IMGT_CDR3:
        return 3

    return None


def renumber_chain(chain: Chain.Chain) -> Chain.Chain:
    """Renumber a chain following IMGT conventions."""
    residues = list(chain.get_residues())

    sequence = ''.join(
        [IUPACData.protein_letters_3to1[res.get_resname().title()] for res in residues if res.id[0] == ' ']
    )

    numbering, chain_type = anarci.number(sequence)

    if not numbering:
        msg = f'Chain ID {chain.id} not identified by ANARCI'
        raise ValueError(msg)

    if len(numbering) == 2:  # noqa: PLR2004
        numbering = numbering[0]
        chain_type = chain_type[0]

        logger.warning(
            ('Multiple possible chain annotations found for chain id %s. ' 'Defaulting to first: %sCHAIN'),
            chain.id,
            chain_type,
        )

    numbering = [(seq_id, insert_code) for (seq_id, insert_code), res_name in numbering if res_name != '-']

    num_residues_not_numbered = len(residues) - len(numbering)

    next_seq_id = numbering[-1][0] + 1

    for _ in range(num_residues_not_numbered):
        numbering.append((next_seq_id, ' '))
        next_seq_id += 1

    output_chain = Chain.Chain(chain.id)
    for (seq_id, insert_code), res in zip(numbering, residues, strict=False):
        new_res = res.copy()
        new_res.id = (new_res.id[0], seq_id, insert_code)
        output_chain.add(new_res)

    return output_chain
