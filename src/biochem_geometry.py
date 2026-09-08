# -*- coding: utf-8 -*-

"""
Biochemical geometry module for Protein-Protein interface analysis.

This module provides the structural cleaning pipeline, 3D spatial distance 
calculations, angle measurements and aromatic ring geometry solvers
(centroids and normal vectors) using Biopython coordinates and NumPy.
"""

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
    """ Filters out water molecules, ligands and heteroatoms.
    
    Args:
        res(Bio.PDB.Residue.Residue): The Biopython Residue object to evaluate.

    Returns:
        bool: True if the residue is a standard amino acid, False if it is water
        a ligand or an active site co-factor.
    """
    res_id = res.get_id()[0]
    return res_id.strip() == ''


def get_clean_atoms(chain):
    """ 
    Extracts heavy atoms from a protein chain ignoring hydrogen atoms and water molecules
    and resolving disordered alternative conformations.
    
    Args: 
        chain(Bio.PDB.Chain.Chain): The Biopython chain object to process.

    Returns:
        list of Bio.PDB.Atom.Atom: A list of clean, non-disordered heavy atoms 
                                   eligible for spatial distance calculations.
    """
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
    """ 
    Calculates the Euclidean distance between two 3D coordinates.
    
    Args:
        coord1 (numpy.ndarray): 3D coordinates of the first atom [x,y,z].
        coord2 (numpy.ndarray): 3D coordinates of the second atom [x,y,z].

    Returns:
    float: The physical Euclidean distance between the 2 points in Angstroms (Å).

    """
    return math.sqrt(
        (coord1[0] - coord2[0])**2 +
        (coord1[1] - coord2[1])**2 +
        (coord1[2] - coord2[2])**2
    )

def calculate_angle(coord_donneur, coord_accepteur, coord):
    """ 
    Calculates the 3D angle in degrees between the antecedent atom,
    the donor atom, and the acceptor atom.

    The mathematical formulation uses the definition of the dot product:
        cos(theta) = (v1 . v2) / (||v1|| * ||v2||)

    Args:
        coord_vertex (numpy.ndarray): Coordinates of the vertex (pivot) atom [x,y,z].
        coord_p1 (numpy.ndarray): Coordinates of the first vector endpoint [x,y,z].
        coord_p2 (numpy.ndarray): Coordinates of the second vector endpoint [x,y,z].

    Returns:
        float: The angle in degrees, ranging from 0.0 to 180.0. Returns
            0.0 if either vector has a magnitude of 0.
    """
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
    """ 
    Extracts the atoms constituting the flat aromatic ring system of a residue.

    Identifies Phe and Tyr (6-membered single rings) or Trp (9-membered indole ring system)
    and returns their corresponding Biopython Atom objects.

    Args:
        residue (Bio.PDB.Residue.Residue): The aromatic residue of interest.

    Returns:
        list of Bio.PDB.Atom.Atom: The Atom objects forming the conjugated plane.
                                   Returns an empty list if the residue is not aromatic.
    """

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
    """ 
    Calculates the 3D centroid and the unit normal vector for an aromatic ring plane.
    
    Args: 
        residue (Bio.PDB.Residue.Residue): The aromatic residue (PHE, TYR or TRP).

    Returns:
        tuple: A tuple that contains:
            - centroid(numpy.ndarray or None): The 3D coordinates [x,y,z] of
            the ring center.
            - normal(numpy.ndarray or None): The normalized unit normal vector [dx,dy,dz]
            perpendicular to the aromatic plane.
            Returns (None, None) if the residue contains fewer than 3 ring atoms.
    """

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