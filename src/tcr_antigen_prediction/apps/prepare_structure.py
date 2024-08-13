'''Generate mesh and compute properties from PDB files.'''
import argparse
import logging
import os

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
    setup_logger(logger, args.log_level)

    base_name = os.path.basename(args.structure).rsplit('.', 1)[0]
    output_name = os.path.join(args.output, base_name)

    if args.chains:
        logger.info('Extracting chains')
        output_name += f"_{''.join(args.chains)}"
        extract_pdb(args.structure, output_name + '.pdb', args.chains)

    logger.info('Re-protonating structure')
    reprotonate(output_name + '.pdb',  output_name + '_protonated.pdb')

    logger.info('Creating surface')
    vertices, faces, _, names, _ = compute_msms(output_name + '_protonated.pdb')

    logger.info('Computing charges')
    vertex_hbond = compute_charges(output_name + '_protonated.pdb', vertices, names)

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
    vertex_charges = compute_apbs(regular_mesh.vertices, output_name + '_protonated.pdb', output_name + '_protonated')

    # TODO decide if this is needed
    # interface = compute_interface(base_name + '_protonated.pdb',
    # ) if args.compute_interface else None

    logging.info('Saving to file')
    save_ply(output_name + '.ply',
             regular_mesh.vertices,
             regular_mesh.faces,
             normals=vertex_normal,
             charges=vertex_charges,
             normalize_charges=True,
             hbond=vertex_hbond,
             hphob=vertex_hphobicity)


if __name__ == '__main__':
    main()
