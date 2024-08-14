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
'''Chemistry constants.'''
import numpy as np

RADII = {
    'N': '1.540000',
    'O': '1.400000',
    'C': '1.740000',
    'H': '1.200000',
    'S': '1.800000',
    'P': '1.800000',
    'Z': '1.39',
    'X': '0.770000',  # Radii of CB or CA in disembodied case.
}
'''Radii for atoms in explicit case.'''

POLAR_HYDROGENS = {
    'ALA': ['H'],
    'GLY': ['H'],
    'SER': ['H', 'HG'],
    'THR': ['H', 'HG1'],
    'LEU': ['H'],
    'ILE': ['H'],
    'VAL': ['H'],
    'ASN': ['H', 'HD21', 'HD22'],
    'GLN': ['H', 'HE21', 'HE22'],
    'ARG': ['H', 'HH11', 'HH12', 'HH21', 'HH22', 'HE'],
    'HIS': ['H', 'HD1', 'HE2'],
    'TRP': ['H', 'HE1'],
    'PHE': ['H'],
    'TYR': ['H', 'HH'],
    'GLU': ['H'],
    'ASP': ['H'],
    'LYS': ['H', 'HZ1', 'HZ2', 'HZ3'],
    'PRO': [],
    'CYS': ['H'],
    'MET': ['H'],
}
'''This  polar hydrogen's names correspond to that of the program Reduce.'''

HBOND_STD_DEV = np.pi / 3

ACCEPTOR_ANGLE_ATOM = {
    'O': 'C',
    'O1': 'C',
    'O2': 'C',
    'OXT': 'C',
    'OT1': 'C',
    'OT2': 'C',
}
'''Dictionary from an acceptor atom to its directly bonded atom on which to compute the angle.'''

ACCEPTOR_PLANE_ATOM = {
    'O': 'CA',
}
'''Dictionary from acceptor atom to a third atom on which to compute the plane.'''

DONOR_ATOM = {
    'H': 'N',
    # Hydrogen bond information.
    # ARG
    # ARG NHX
    # Angle: NH1, HH1X, point and NH2, HH2X, point 180 degrees.
    # radii from HH: radii[H]
    # ARG NE
    # Angle: ~ 120 NE, HE, point, 180 degrees
    'HH11': 'NH1',
    'HH12': 'NH1',
    'HH21': 'NH2',
    'HH22': 'NH2',
    'HE': 'NE',
}
'''Dictionary from an H atom to its donor atom.'''

# ASN
# Angle ND2,HD2X: 180
# Plane: CG,ND2,OD1
# Angle CG-OD1-X: 120
DONOR_ATOM['HD21'] = 'ND2'
DONOR_ATOM['HD22'] = 'ND2'
# ASN Acceptor
ACCEPTOR_ANGLE_ATOM['OD1'] = 'CG'
ACCEPTOR_PLANE_ATOM['OD1'] = 'CB'

# ASP
# Plane: CB-CG-OD1
# Angle CG-ODX-point: 120
ACCEPTOR_ANGLE_ATOM['OD2'] = 'CG'
ACCEPTOR_PLANE_ATOM['OD2'] = 'CB'

# GLU
# PLANE: CD-OE1-OE2
# ANGLE: CD-OEX: 120
# GLN
# PLANE: CD-OE1-NE2
# Angle NE2,HE2X: 180
# ANGLE: CD-OE1: 120
DONOR_ATOM['HE21'] = 'NE2'
DONOR_ATOM['HE22'] = 'NE2'
ACCEPTOR_ANGLE_ATOM['OE1'] = 'CD'
ACCEPTOR_ANGLE_ATOM['OE2'] = 'CD'
ACCEPTOR_PLANE_ATOM['OE1'] = 'CG'
ACCEPTOR_PLANE_ATOM['OE2'] = 'CG'

# HIS Acceptors: ND1, NE2
# Plane ND1-CE1-NE2
# Angle: ND1-CE1 : 125.5
# Angle: NE2-CE1 : 125.5
ACCEPTOR_ANGLE_ATOM['ND1'] = 'CE1'
ACCEPTOR_ANGLE_ATOM['NE2'] = 'CE1'
ACCEPTOR_PLANE_ATOM['ND1'] = 'NE2'
ACCEPTOR_PLANE_ATOM['NE2'] = 'ND1'

# HIS Donors: ND1, NE2
# Angle ND1-HD1 : 180
# Angle NE2-HE2 : 180
DONOR_ATOM['HD1'] = 'ND1'
DONOR_ATOM['HE2'] = 'NE2'

# TRP Donor: NE1-HE1
# Angle NE1-HE1 : 180
DONOR_ATOM['HE1'] = 'NE1'

# LYS Donor NZ-HZX
# Angle NZ-HZX : 180
DONOR_ATOM['HZ1'] = 'NZ'
DONOR_ATOM['HZ2'] = 'NZ'
DONOR_ATOM['HZ3'] = 'NZ'

# TYR acceptor OH
# Plane: CE1-CZ-OH
# Angle: CZ-OH 120
ACCEPTOR_ANGLE_ATOM['OH'] = 'CZ'
ACCEPTOR_PLANE_ATOM['OH'] = 'CE1'

# TYR donor: OH-HH
# Angle: OH-HH 180
DONOR_ATOM['HH'] = 'OH'
ACCEPTOR_PLANE_ATOM['OH'] = 'CE1'

# SER acceptor:
# Angle CB-OG-X: 120
ACCEPTOR_ANGLE_ATOM['OG'] = 'CB'

# SER donor:
# Angle: OG-HG-X: 180
DONOR_ATOM['HG'] = 'OG'

# THR acceptor:
# Angle: CB-OG1-X: 120
ACCEPTOR_ANGLE_ATOM['OG1'] = 'CB'

# THR donor:
# Angle: OG1-HG1-X: 180
DONOR_ATOM['HG1'] = 'OG1'
