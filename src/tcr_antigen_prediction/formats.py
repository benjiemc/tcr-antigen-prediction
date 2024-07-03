from Bio.PDB import PDBParser

from tcr_antigen_prediction.chemistry import radii, polarHydrogens


def output_pdb_as_xyzrn(pdb_filename, xyzrn_filename):
    """
        pdb_filename: input pdb filename
        xyzrn_filename: output in xyzrn format.
    """
    parser = PDBParser()
    struct = parser.get_structure(pdb_filename, pdb_filename)
    with open(xyzrn_filename, "w") as outfile:
        for atom in struct.get_atoms():
            name = atom.get_name()
            residue = atom.get_parent()
            # Ignore hetatms.
            if residue.get_id()[0] != " ":
                continue

            resname = residue.get_resname()
            chain = residue.get_parent().get_id()
            atomtype = name[0]

            color = "Green"
            coords = None
            if atomtype in radii and resname in polarHydrogens:
                if atomtype == "O":
                    color = "Red"
                if atomtype == "N":
                    color = "Blue"
                if atomtype == "H":
                    if name in polarHydrogens[resname]:
                        color = "Blue"  # Polar hydrogens
                coords = "{:.06f} {:.06f} {:.06f}".format(
                    atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
                )
                insertion = "x"
                if residue.get_id()[2] != " ":
                    insertion = residue.get_id()[2]
                full_id = "{}_{:d}_{}_{}_{}_{}".format(
                    chain, residue.get_id()[1], insertion, resname, name, color
                )
            if coords is not None:
                outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")
