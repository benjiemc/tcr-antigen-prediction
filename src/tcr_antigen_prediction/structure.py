# This file is part of the MaSIF project.
#
# Copyright 2019 - Gainza P, Sverrisson F, Monti F, Rodola, Bronstein MM, Correia BE
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
'''Functions and classes for interacting with PDB structures.'''
from subprocess import Popen, PIPE
from typing import Optional, Set, Iterable

import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Selection, StructureBuilder, Select, Structure, Model
from Bio.SeqUtils import IUPACData
PROTEIN_LETTERS = [x.upper() for x in IUPACData.protein_letters_3to1.keys()]


# Exclude disordered atoms.
class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == 'A' or atom.get_altloc() == '1'


def find_modified_amino_acids(path):
    '''Contributed by github user jomimc - find modified amino acids in the PDB (e.g. MSE)'''
    res_set = set()

    for line in open(path, 'r'):
        if line[:6] == 'SEQRES':
            for res in line.split()[4:]:
                res_set.add(res)

    for res in list(res_set):
        if res in PROTEIN_LETTERS:
            res_set.remove(res)

    return res_set


def extract_pdb(infilename, outfilename, chain_ids=None):
    # extract the chain_ids from infilename and save in outfilename.
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(infilename, infilename)
    model = Selection.unfold_entities(struct, 'M')[0]

    # Select residues to extract and build new structure
    struct_builder = StructureBuilder.StructureBuilder()

    struct_builder.init_structure('output')
    struct_builder.init_seg(' ')
    struct_builder.init_model(0)

    output_structure = struct_builder.get_structure()

    # Load a list of non-standard amino acid names -- these are typically listed under HETATM, so they would be
    # typically ignored by the orginal algorithm
    modified_amino_acids = find_modified_amino_acids(infilename)

    for chain in model:
        if chain_ids is None or chain.get_id() in chain_ids:
            struct_builder.init_chain(chain.get_id())

            for residue in chain:
                het = residue.get_id()

                if het[0] == ' ':
                    output_structure[0][chain.get_id()].add(residue)

                elif het[0][-3:] in modified_amino_acids:
                    output_structure[0][chain.get_id()].add(residue)

    # Output the selected residues
    pdbio = PDBIO()
    pdbio.set_structure(output_structure)
    pdbio.save(outfilename, select=NotDisordered())


def reprotonate(path: str, out_path: str):
    '''Remove hydrogens (if any) and re-protonate a structure.'''
    # Remove hydrogens
    args = ['reduce', '-Trim', path]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, _ = p2.communicate()

    with open(out_path, 'w') as outfile:
        outfile.write(stdout.decode('utf-8').rstrip())

    # Re-add hydrogens
    args = ['reduce', '-HIS', out_path]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, _ = p2.communicate()

    with open(out_path, 'w') as outfile:
        outfile.write(stdout.decode('utf-8'))


def get_sequence(structure: Structure.Structure, chain_id: str, residue_range: Optional[Set[int]] = None) -> str:
    '''Get the sequence of amino acids from a biopython strcture.

    Args:
        structure: Biopython structure.
        chain_id: ID of the chain to get the sequence from.
        residue_range: Range of residues to include in the sequence (Optional). Default is to include them all.

    Returns:
        Amino acid sequence as one-letter codes.

    '''
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
    '''Get the header lines from the contents of a pdb file.'''
    header = []
    for line in pdb_contents.split('\n'):
        record_type = line[0:6]
        if record_type in ('ATOM  ', 'HETATM', 'MODEL '):
            break

        header.append(line)

    return '\n'.join(header)


def extract_chains(structure: Structure.Structure, chains: Iterable[str]) -> Structure.Structure:
    '''Get only the selected chain from a PDB structure.'''
    new_structure = Structure.Structure(structure.id)

    for model in structure:
        new_model = Model.Model(model.id)
        for chain in model:
            if chain.id in chains:
                new_model.add(chain.copy())

        new_structure.add(new_model)

    return new_structure


def bio_to_pandas(structure: Structure.Structure) -> pd.DataFrame:
    '''Convert a biopython structure to a pandas dataframe.'''
    records = []

    multiple_models = len(structure) > 1

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    records.append({
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
                    })

                    if multiple_models:
                        records[-1]['model_index'] = model.id

        column_names = ['record_type', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
                        'residue_seq_id', 'residue_insert_code', 'pos_x', 'pos_y', 'pos_z', 'occupancy',
                        'b_factor', 'element', 'charge']
    if multiple_models:
        column_names.append('model_index')

    return pd.DataFrame(records, columns=column_names)
