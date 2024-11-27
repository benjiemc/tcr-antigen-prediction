"""Functions for adding missing residue information to sequences."""

import logging
import re
from collections.abc import Iterable

from Bio.PDB import Structure
from Bio.SeqUtils import IUPACData

from tcr_antigen_prediction.aligners import align_sequences
from tcr_antigen_prediction.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, IMGT_VARIABLE_DOMAIN

logger: logging.Logger = logging.getLogger(__name__)


def get_missing_residues(header: str) -> list[dict]:
    """Extract all missing residues from the header of a pdb file.

    Returns residue name, chain id, residue/sequence id, and insert id if available.
    """
    lines = re.findall(r'^REMARK 465\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?', header, flags=re.MULTILINE)
    return [
        {'residue_name': res_name, 'chain_id': chain, 'residue_seq_id': int(seq_id), 'residue_insert_code': insert_id}
        for res_name, chain, seq_id, insert_id in lines
    ]


def get_missing_atoms(header: str) -> list[dict]:
    """Extract residues with missing atoms from the header of a pdb file."""
    lines = re.findall(r'^REMARK 470\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?', header, flags=re.MULTILINE)
    return [
        {'residue_name': res_name, 'chain_id': chain, 'residue_seq_id': int(seq_id), 'residue_insert_code': insert_id}
        for res_name, chain, seq_id, insert_id in lines
    ]


def compare_sequences_for_missing(
    sequence: list[list], annotated_raw_sequence: list[list], numbering: set[int]
) -> bool:
    """Check if there are missing residues in the numbered regions."""
    seq = ''.join([res[0] for res in sequence])
    raw_seq = ''.join([res[0] for res in annotated_raw_sequence])

    alignment, _ = align_sequences(raw_seq, seq)

    raw_index = 0
    index = 0
    current_seq_id = 0

    for raw_res, res in alignment:
        if raw_res == '-':
            index += 1
            continue

        if res != '-':
            current_seq_id = sequence[index][1]

        if annotated_raw_sequence[raw_index][3] == 'missing' and current_seq_id in numbering:
            return False

        if res == '-':
            raw_index += 1
            continue

        index += 1
        raw_index += 1

    return True


def trim_start_and_end(annotated_raw_sequence: list[list]) -> list[list]:
    """Remove missing resiudes at start and end of the sequence."""
    start_index = 0
    for res in annotated_raw_sequence:
        if res[3] != 'missing':
            break

        start_index += 1

    end_index = len(annotated_raw_sequence)
    for res in reversed(annotated_raw_sequence):
        if res[3] != 'missing':
            break

        end_index -= 1

    return annotated_raw_sequence[start_index:end_index]


def screen_chain(
    structure: Structure.Structure,
    raw_structure: Structure.Structure,
    chain_id: str,
    missing_residues: list[dict],
    missing_atoms: list[dict],
    numbering: set[int],
) -> bool:
    """Screen chain for missing residues or atoms."""
    missing_residues_on_chain = [
        [
            IUPACData.protein_letters_3to1[entity['residue_name'].title()],
            entity['residue_seq_id'],
            entity['residue_insert_code'],
            'missing',
        ]
        for entity in missing_residues
        if entity['chain_id'] == chain_id
    ]

    missing_atoms_on_chain = [
        [
            IUPACData.protein_letters_3to1[entity['residue_name'].title()],
            entity['residue_seq_id'],
            entity['residue_insert_code'],
        ]
        for entity in missing_atoms
        if entity['chain_id'] == chain_id
    ]

    raw_sequence_info = [
        [IUPACData.protein_letters_3to1[res.get_resname().title()], res.id[1], res.id[2].strip(), 'not-missing']
        for res in raw_structure[0][chain_id].get_residues()
        if res.id[0] == ' '
    ]
    imgt_sequence_info = [
        [IUPACData.protein_letters_3to1[res.get_resname().title()], res.id[1], res.id[2].rstrip()]
        for res in structure[0][chain_id].get_residues()
        if res.id[0] == ' '
    ]

    combined_raw_residues = sorted(missing_residues_on_chain + raw_sequence_info, key=lambda res: (res[1], res[2]))
    combined_raw_residues = trim_start_and_end(combined_raw_residues)

    for res in combined_raw_residues:
        if res[:3] in missing_atoms_on_chain:
            res[3] = 'missing'

    return compare_sequences_for_missing(imgt_sequence_info, combined_raw_residues, numbering)


def screen_tcr_variable_domain(
    structure: Structure.Structure,
    raw_structure: Structure.Structure,
    tcr_chain_ids: Iterable[str],
    missing_residues: list[dict],
    missing_atoms: list[dict],
) -> bool:
    """Check if there are missing residues or atoms in the TCR variable domain.

    Residues are allowed to be missing at the start or end of the structure.
    TODO: allow multiple models

    Args:
        structure: imgt numbered biopython structure
        raw_structure: pdb biopython structure
        tcr_chain_ids: chain ids for the tcr alpha and beta chains
        missing_residues: list of missing residues in the structure formatted as dicts with the following keys
                          'residue_name', 'chain_id', 'residue_seq_id', 'residue_insert_code'
        missing_atoms: list of residues missing atoms in the structure formatted as dicts with the following keys
                          'residue_name', 'chain_id', 'residue_seq_id', 'residue_insert_code'

    Returns:
        True is there are no missing atoms or residues in the variable domain and False if there are

    """
    for chain_id in tcr_chain_ids:
        if not screen_chain(structure, raw_structure, chain_id, missing_residues, missing_atoms, IMGT_VARIABLE_DOMAIN):
            return False

    return True


def screen_pmhc_abd(
    structure: Structure.Structure,
    raw_structure: Structure.Structure,
    antigen_chain_id: str,
    mhc_chains: Iterable[str],
    mhc_type: str,
    missing_residues: list[dict],
    missing_atoms: list[dict],
) -> bool:
    """Screen the pMHC antigen binding domain for missing atoms or residues.

    This includes the MHC antigen binding domain and the antigen. Residues are allowed to be missing at the start or end
    of the MHC antigen binding domain.

    TODO: allow multiple models

    Args:
        structure: imgt numbered biopython structure
        raw_structure: pdb biopython structure
        antigen_chain_id: id for the antigen
        mhc_chains: chain ids for the MHC domain
        mhc_type: type of mhc molecule, can either be MH1 or MH2
        missing_residues: list of missing residues in the structure formatted as dicts with the following keys
                          'residue_name', 'chain_id', 'residue_seq_id', 'residue_insert_code'
        missing_atoms: list of residues missing atoms in the structure formatted as dicts with the following keys
                          'residue_name', 'chain_id', 'residue_seq_id', 'residue_insert_code'

    Returns:
        True if there are not missing residues or atoms, otherwise False

    Raises:
        ValueError if the mhc type is not 'MH1' or 'MH2'

    """
    if antigen_chain_id in ({res['chain_id'] for res in missing_residues} | {res['chain_id'] for res in missing_atoms}):
        return False

    if mhc_type == 'MH1':
        if not screen_chain(structure, raw_structure, mhc_chains[0], missing_residues, missing_atoms, IMGT_MH1_ABD):
            return False

    elif mhc_type == 'MH2':
        for chain_id in mhc_chains:
            if not screen_chain(structure, raw_structure, chain_id, missing_residues, missing_atoms, IMGT_MH2_ABD):
                return False

    else:
        msg = 'Invalid MHC type, must be MH1 or MH2'
        raise ValueError(msg)

    return True
