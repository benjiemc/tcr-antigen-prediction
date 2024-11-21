"""Functions and classes for interacting with PDB structures."""

from collections.abc import Iterable

import pandas as pd
from Bio.PDB import Chain, Model, Structure
from Bio.SeqUtils import IUPACData

PROTEIN_LETTERS: list[str] = [x.upper() for x in IUPACData.protein_letters_3to1]
'''Amino acid one letter codes.'''


def get_sequence(structure: Structure.Structure, chain_id: str, residue_range: set[int] | None = None) -> str:
    """Get the sequence of amino acids from a biopython strcture.

    Args:
        structure: Biopython structure.
        chain_id: ID of the chain to get the sequence from.
        residue_range: Range of residues to include in the sequence (Optional). Default is to include them all.

    Returns:
        Amino acid sequence as one-letter codes.

    """
    chain = structure[0][chain_id]
    sequence = []

    for residue in chain:
        if residue.id[0] != ' ':
            continue

        if residue_range is not None and residue.id[1] not in residue_range:
            continue

        sequence.append(residue.get_resname())

    sequence = [IUPACData.protein_letters_3to1[res.title()] for res in sequence]

    return ''.join(sequence)


def get_header(pdb_contents: str) -> str:
    """Get the header lines from the contents of a pdb file."""
    header = []
    for line in pdb_contents.split('\n'):
        record_type = line[0:6]
        if record_type in ('ATOM  ', 'HETATM', 'MODEL '):
            break

        header.append(line)

    return '\n'.join(header)


def extract_chains(structure: Structure.Structure, chains: Iterable[str]) -> Structure.Structure:
    """Get only the selected chain from a PDB structure."""
    new_structure = Structure.Structure(structure.id)

    for model in structure:
        new_model = Model.Model(model.id)
        for chain in model:
            if chain.id in chains:
                new_model.add(chain.copy())

        new_structure.add(new_model)

    return new_structure


def bio_to_pandas(structure: Structure.Structure) -> pd.DataFrame:
    """Convert a biopython structure to a pandas dataframe."""
    records = []

    multiple_models = len(structure) > 1

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    records.append(
                        {
                            'record_type': 'ATOM' if residue.id[0] == ' ' else 'HETATM',
                            'atom_number': atom.serial_number,
                            'atom_name': atom.id,
                            'alt_loc': atom.altloc if atom.altloc != ' ' else None,
                            'residue_name': residue.resname,
                            'chain_id': chain.id,
                            'residue_seq_id': residue.id[1],
                            'residue_insert_code': residue.id[2] if residue.id[2] != ' ' else None,
                            'pos_x': atom.coord[0],
                            'pos_y': atom.coord[1],
                            'pos_z': atom.coord[2],
                            'occupancy': atom.occupancy,
                            'b_factor': atom.bfactor,
                            'element': atom.element,
                            'charge': atom.get_charge(),
                        }
                    )

                    if multiple_models:
                        records[-1]['model_index'] = model.id

        column_names = [
            'record_type',
            'atom_number',
            'atom_name',
            'alt_loc',
            'residue_name',
            'chain_id',
            'residue_seq_id',
            'residue_insert_code',
            'pos_x',
            'pos_y',
            'pos_z',
            'occupancy',
            'b_factor',
            'element',
            'charge',
        ]
    if multiple_models:
        column_names.append('model_index')

    return pd.DataFrame(records, columns=column_names)


def replace_chain(structure: Structure.Structure, new_chain: Chain.Chain) -> Structure.Structure:
    """Replace chain with a new chain in a PDB structure (does not modify original structure)."""
    output_structure = structure.copy()
    chain_id = new_chain.id

    for model in output_structure:
        model.detach_child(chain_id)
        model.add(new_chain)

    return output_structure
