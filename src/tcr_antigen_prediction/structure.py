'''Functions and classes for interacting with PDB structures.'''
from subprocess import Popen, PIPE

from Bio.PDB import PDBParser, PDBIO, Selection, StructureBuilder, Select
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
