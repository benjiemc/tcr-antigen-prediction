import argparse
import logging
import os

import pymesh

from tcr_antigen_prediction.protonate import reprotonate
from tcr_antigen_prediction.structure import extractPDB
from tcr_antigen_prediction.triangulate import (fix_mesh,
                                                computeMSMS,
                                                computeCharges,
                                                computeHydrophobicity,
                                                assignChargesToNewMesh,
                                                compute_normal,
                                                computeAPBS,
                                                save_ply)

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(prog='PrepareStructure')

parser.add_argument('structure', help='Path to pdb structure.')
parser.add_argument('--output', '-o', help='Path to output processed structure')
parser.add_argument('--chains', default=None, nargs='+', help='Chains to use from protein.')
parser.add_argument('--compute-interface', action='store_true')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    args = parser.parse_args()

    base_name = os.path.basename(args.structure).rsplit('.', 1)[0]
    output_name = os.path.join(args.output, base_name)

    if args.chains:
        logger.info('Extracting chains')
        output_name += f"_{''.join(args.chains)}"
        extractPDB(args.structure, output_name + '.pdb', args.chains)

    logger.info('Re-protonating structure')
    reprotonate(output_name + '.pdb',  output_name + '_protonated.pdb')

    logger.info('Creating surface')
    vertices, faces, normals, names, areas = computeMSMS(output_name + '_protonated.pdb')

    logger.info('Computing charges')
    vertex_hbond = computeCharges(output_name + '_protonated.pdb', vertices, names)
    logger.info('Computing hydrophobicity')
    vertex_hphobicity = computeHydrophobicity(names)

    logger.info('Creating mesh')
    mesh = pymesh.form_mesh(vertices, faces)
    logger.info('Regularizing mesh')
    regular_mesh = fix_mesh(mesh)

    logger.info('Computing normals')
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices, vertex_hbond)
    vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices, vertex_hphobicity)
    vertex_charges = computeAPBS(regular_mesh.vertices, output_name + '_protonated.pdb', output_name + '_protonated')

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
