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
'''Generate mesh and compute properties from PDB files.'''
import argparse
import logging
import os
import shutil
import tempfile

import pymesh

from tcr_antigen_prediction.apps._log import add_logging_arguments, setup_logger
from tcr_antigen_prediction.io import save_ply
from tcr_antigen_prediction.structure import extract_pdb, reprotonate
from tcr_antigen_prediction.triangulate import (fix_mesh,
                                                compute_msms,
                                                compute_charges,
                                                compute_hydrophobicity,
                                                assign_charges_to_new_mesh,
                                                compute_normal,
                                                compute_apbs)

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog='PrepareStructure')

parser.add_argument('structure', help='Path to pdb structure.')
parser.add_argument('--output', '-o', help='Path to output processed structure')
parser.add_argument('--chains', default=None, nargs='+', help='Chains to use from protein.')
parser.add_argument('--compute-interface', action='store_true')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    base_name = os.path.basename(args.structure).rsplit('.', 1)[0]

    with tempfile.TemporaryDirectory() as temp_dir:

        if args.chains:
            logger.info('Extracting chains')
            base_name += f"_{''.join(args.chains)}"
            extract_pdb(args.structure, os.path.join(temp_dir, base_name + '.pdb'), args.chains)

        else:
            logger.debug('Copying PDB file to workspace')
            shutil.copy(args.structure, os.path.join(temp_dir, base_name + '.pdb'))

        logger.info('Re-protonating structure')
        reprotonate(os.path.join(temp_dir, base_name + '.pdb'),
                    os.path.join(temp_dir, base_name + '_protonated.pdb'))

        logger.info('Creating surface')
        vertices, faces, _, names, _ = compute_msms(os.path.join(temp_dir, base_name + '_protonated.pdb'))

        logger.info('Computing charges')
        vertex_hbond = compute_charges(os.path.join(temp_dir, base_name + '_protonated.pdb'), vertices, names)

        logger.info('Computing hydrophobicity')
        vertex_hphobicity = compute_hydrophobicity(names)

        logger.info('Creating mesh')
        mesh = pymesh.form_mesh(vertices, faces)

        logger.debug('Regularizing mesh')
        regular_mesh = fix_mesh(mesh)

        logger.info('Computing normals')
        vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

        vertex_hbond = assign_charges_to_new_mesh(regular_mesh.vertices, vertices, vertex_hbond)
        vertex_hphobicity = assign_charges_to_new_mesh(regular_mesh.vertices, vertices, vertex_hphobicity)
        vertex_charges = compute_apbs(regular_mesh.vertices, os.path.join(temp_dir, base_name + '_protonated.pdb'))

    # TODO decide if this is needed
    # interface = compute_interface(base_name + '_protonated.pdb',
    # ) if args.compute_interface else None

    logging.info('Saving to file')
    save_ply(os.path.join(args.output, base_name) + '.ply',
             regular_mesh.vertices,
             regular_mesh.faces,
             normals=vertex_normal,
             charges=vertex_charges,
             normalize_charges=True,
             hbond=vertex_hbond,
             hphob=vertex_hphobicity)


if __name__ == '__main__':
    main()
