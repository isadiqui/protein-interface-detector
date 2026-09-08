# -*- coding: utf-8 -*-


import math
import numpy as np

# Global dictionaries for biochemical classification
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "MET", "PRO"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
RES_POSITIFS = {"LYS", "ARG", "HIS"}
RES_NEGATIFS = {"ASP", "GLU"}
N_CATION = {"NZ", "NH1", "NH2", "ND1", "NE2"}
O_ANION = {"OD1", "OD2", "OE1", "OE2"}
AROMATIC_RES = {"PHE", "TYR", "TRP"}

def is_res(res):
    """ Filters out water molecules and heteroatoms. """
    res_id = res.get_id()[0]
    return res_id.strip() == ''

def get_clean_atoms(chain):
    """ Extracts heavy atoms from a chain ignoring hydrogen atoms and water molecules
    and resolving disordered alternative conformations. """
    clean_atoms = []
    for res in chain:
        if not is_res(res):
            continue
        
        for atom in res:
            # Ignore hydrogen atoms
            if atom.element == 'H':
                continue

            # Resolve disordered alternative conformations
            if atom.is_disordered():
                clean_atom = atom.selected_child
            else:
                clean_atom = atom
            
            clean_atoms.append(clean_atom)
        
    return clean_atoms

def calculate_distance(coord1, coord2):
    """ Calculates the Euclidean distance between two 3D coordinates. """
    return math.sqrt(
        (coord1[0] - coord2[0])**2 +
        (coord1[1] - coord2[1])**2 +
        (coord1[2] - coord2[2])**2
    )

def calculate_angle(coord_donneur, coord_accepteur, coord):
    """ Calculates the angle in degrees between the antecedent atom,
    the donor atom, and the acceptor atom. """
    v1 = coord - coord_donneur
    v2 = coord_accepteur - coord_donneur
    dot_product = sum(v1 * v2)
    norm_v1 = math.sqrt(sum(v1**2))
    norm_v2 = math.sqrt(sum(v2**2))

    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0
    cos_angle = max(-1.0, min(1.0, dot_product / (norm_v1 * norm_v2)))
    return math.degrees(math.acos(cos_angle))


def get_arom_atoms(residue):
    """ Returns the list of atoms forming the aromatic ring of the residue. """
    resname = residue.get_resname()
    if resname in {"PHE", "TYR"}:
        # 6 membered phenyl/phenol ring
        ring_atom_names = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
    elif resname == "TRP":
        # 9 membered double ring system
        ring_atom_names = {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
    else:
        return []

    return [residue[name] for name in ring_atom_names if name in residue]


def calculate_centroid_normal(residue):
    """ Calculates the 3D centroid and the unit normal vector for an aromatic ring plane. """
    atoms = get_arom_atoms(residue)
    if len(atoms) < 3:
        return None, None

    coords = np.array([atom.get_coord() for atom in atoms])
    centroid = np.mean(coords, axis=0)

    # Select 3 tripod atoms separated by roughly 1/3 of the ring length
    p0, p1, p2 = coords[0], coords[len(coords)//3], coords[2*len(coords)//3]
    v1 = p1 - p0
    v2 = p2 - p0

    # Vector cross-product to get normal vector, then normalize it to 1.0 length
    normal = np.cross(v1, v2)
    normal = normal / np.linalg.norm(normal)

    return centroid, normal