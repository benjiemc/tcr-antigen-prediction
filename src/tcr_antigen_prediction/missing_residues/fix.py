"""Functions and classes for modelling missing residues in structures using MODELLER.

Requires modeller to be available can be downloaded here: https://salilab.org/modeller/

"""

import contextlib
import io
import logging
import os
import sys
from tempfile import TemporaryDirectory

from Bio.PDB import PDBIO, Chain, PDBParser, Structure

logger = logging.getLogger(__name__)

try:
    from modeller import Environ, Selection, log
    from modeller.automodel import AutoModel

except ImportError:
    logger.exception('Modeller not found. Please follow install instructions here: https://salilab.org/modeller/')
    sys.exit(1)

_ALIGNMENT_TEMPLATE = '''>P1;target
sequence:target:::::::0.00:0.00
{sequence}*

>P1;struct
structure:struct:{start}:{chain_id}:::::0.00:0.00
{template}*'''


class MissingResiduesModel(AutoModel):
    """AutoModel class for modelling missing resiudes.

    Only the missing residues can be optimized.

    Attributes:
        missing_residue_ranges: list of ranges (inclusive) of residues missing in the template.

    """

    def __init__(self, missing_residue_ranges: list[tuple[int, int]], *args: tuple, **kwargs: dict) -> None:
        """Initialize the class.

        Args:
            missing_residue_ranges: list of ranges (inclusive) of residues missing in the template.
            *args: arguments for AutoModel class
            **kwargs: key-worded arguments for AutoModel class

        """
        super().__init__(*args, **kwargs)
        self.missing_residue_ranges = missing_residue_ranges

    def select_atoms(self) -> Selection:
        """Atoms selected for optimization are only in the missing residues range."""
        return Selection(*[self.residue_range(f'{start}:A', f'{end}:A') for start, end in self.missing_residue_ranges])


def get_ranges(template: str) -> list[tuple[int, int]]:
    """Get missing residue ranges in template sequence.

    Example:
        full sequence: SRPWFLEYCPTKTQPLEHHNLLVCSVSDFYPGNIEVRWFRNGKEEKTGIVSTGLVRNGDWTFQTLVMLETVPQSGEVYTCQVEHPSLTDPVTV
        template:      --PWFLEYCPT----L----L-VCSVSDFYPGNIEVRWF---KEEKTGIVSTGLVRNGDWTFQTLVMLET-------YTCQVEHPSLTDPVTV
                       **         **** **** *                 ***                            *******

        template = '--PWFLEYCPT----L----L-VCSVSDFYPGNIEVRWF---KEEKTGIVSTGLVRNGDWTFQTLVMLET-------YTCQVEHPSLTDPVTV'
        get_ranges(template)
        [(1, 2), (11, 14), (16, 19), (21, 21), (38, 40), (68, 74)]

    Args:
        template: sequence with dashes representing missing residues

    Returns:
        ranges (indexed from 1) with inclusive stop and start values

    """
    ranges = []

    range_start = None
    range_end = None

    for pos, amino_acid in enumerate(template, 1):
        if amino_acid == '-' and range_start is None:
            range_start = pos

        elif amino_acid != '-' and range_start is not None:
            range_end = pos - 1

            ranges.append((range_start, range_end))

            range_start = None
            range_end = None

    if range_start is not None:
        ranges.append((range_start, len(template)))

    return ranges


def predict_missing_residues(
    alignment: list[tuple[str, str]], structure: Structure.Structure, chain_id: str
) -> Chain.Chain:
    """Use MODELLER to fill in missing residues on a protein chain.

    Args:
        alignment: aligned complete sequence (including missing residues) and template sequence from structure
        structure: template PDB structure
        chain_id: ID of the chain to fill in missing residues

    Returns:
        new chain with filled in residues

    """
    sequence = ''.join([res for res, _ in alignment])
    template = ''.join([res for _, res in alignment])

    missing_residue_ranges = get_ranges(template)

    residues = [res for res in structure[0][chain_id].get_residues() if res.id[0] == ' ']
    start_pos = f'{residues[0].id[1]}{residues[0].id[2].strip()}'

    with TemporaryDirectory() as temp_dir:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(temp_dir, 'struct.pdb'))

        align_path = os.path.join(temp_dir, 'alignment.ali')
        with open(align_path, 'w') as fh:
            fh.write(
                _ALIGNMENT_TEMPLATE.format(sequence=sequence, template=template, chain_id=chain_id, start=start_pos)
            )

        with contextlib.chdir(temp_dir):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                log.none()
                env = Environ()
                env.io.atom_files_directory = ['.']

                prediction = MissingResiduesModel(
                    missing_residue_ranges,
                    env,
                    alnfile='alignment.ali',
                    knowns='struct',
                    sequence='target',
                    root_name='output_struct',
                )

                prediction.starting_model = 1
                prediction.ending_model = 1

                prediction.make()

            logger.debug(output.getvalue())

        pdb_parser = PDBParser(QUIET=True)
        fixed_structure = pdb_parser.get_structure('', os.path.join(temp_dir, 'output_struct.B99990001.pdb'))

    fixed_chain = fixed_structure[0]['A']
    fixed_chain.detach_parent()

    fixed_chain.id = chain_id

    return fixed_chain
