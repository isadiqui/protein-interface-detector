#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
Protein Interface Contact Analyzer (PICA)
Author: Master Bioinformatique, Université de Paris Cité
Project: UE Programmation et Gestion de Projet - Short Project
================================================================================
This script parses a PDB file, cleans the structure, and identifies/characterizes
four types of non-covalent interactions (hydrophobic, ionic/salt bridges, 
hydrogen bonds, and aromatic pi-pi stacking) at the interface of two protein chains.
It generates a structured report, exports results to CSV, and automatically writes
a PyMOL visualization script (.pml) for 3D analysis.
================================================================================
"""

import os
import sys
import math
import csv
import argparse
import numpy as np
from Bio.PDB import PDBParser

# --- GLOBAL DICTIONARIES FOR BIOCHEMICAL CLASSIFICATION ---
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "MET", "PRO"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

RES_POSITIFS = {"LYS", "ARG", "HIS"}
RES_NEGATIFS = {"ASP", "GLU"}
N_CATION = {"NZ", "NH1", "NH2", "ND1", "NE2"}
O_ANION = {"OD1", "OD2", "OE1", "OE2"}

AROMATIC_RES = {"PHE", "TYR", "TRP"}


def is_res(res):
    """
    Filters out water molecules and heteroatoms.
    In Biopython, a residue ID starting with ' ' (blank) indicates a standard
    amino acid residue, whereas 'W' indicates water and 'H_' represents ligands/heteroatoms.
    """
    res_id = res.get_id()[0]
    return res_id.strip() == ''


def get_clean_atoms(chain):
    """
    Extracts heavy atoms from a chain, ignoring hydrogen atoms and water molecules,
    and resolving alternative disordered conformations to avoid coordinates conflicts.
    """
    clean_atoms = []
    for res in chain:
        if not is_res(res):
            continue
        for atom in res:
            # Ignore hydrogens
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
    """Calculates the Euclidean distance between two 3D coordinates."""
    return math.sqrt(
        (coord1[0] - coord2[0])**2 + 
        (coord1[1] - coord2[1])**2 + 
        (coord1[2] - coord2[2])**2
    )


def calculate_angle(coord_donneur, coord_accepteur, coord):
    """
    Calculates the angle in degrees between the antecedent atom, 
    the donor atom, and the acceptor atom.
    """
    v1 = coord - coord_donneur
    v2 = coord_accepteur - coord_donneur
    dot_product = sum(v1 * v2)
    norm_v1 = math.sqrt(sum(v1 ** 2))
    norm_v2 = math.sqrt(sum(v2 ** 2))
    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0
    cos_angle = max(-1.0, min(1.0, dot_product / (norm_v1 * norm_v2)))
    return math.degrees(math.acos(cos_angle))


def get_aromatic_atoms(residue):
    """Returns the list of atoms forming the flat aromatic ring of the residue."""
    resname = residue.get_resname()
    if resname in {"PHE", "TYR"}:
        # 6-membered phenyl/phenol ring
        ring_atom_names = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
    elif resname == "TRP":
        # 9-membered double indole ring system
        ring_atom_names = {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
    else:
        return []
    return [residue[name] for name in ring_atom_names if name in residue]


def calculate_centroid_normal(residue):
    """
    Calculates the 3D centroid and the unit normal vector for an aromatic ring plane.
    Uses three spread-out ring atoms to stably define the plane's orientation.
    """
    atoms = get_aromatic_atoms(residue)
    if len(atoms) < 3:
        return None, None
        
    coords = np.array([atom.get_coord() for atom in atoms])
    centroid = np.mean(coords, axis=0)
    
    # Select three tripod atoms separated by roughly 1/3 of the ring length
    p0 = coords[0]
    p1 = coords[len(coords) // 3]
    p2 = coords[2 * len(coords) // 3]
    
    v1 = p1 - p0
    v2 = p2 - p0
    
    # Vector cross-product to get normal vector, then normalize it to 1.0 length
    normal = np.cross(v1, v2)
    norm_val = np.linalg.norm(normal)
    if norm_val == 0:
        return None, None
    normal = normal / norm_val
    
    return centroid, normal


def write_pymol_script(filepath, pdb_path, chain_A_id, chain_B_id, matches):
    """
    Generates a PyMOL script (.pml) to visually highlight the identified 
    interface residues and depict specific interactions.
    """
    with open(filepath, 'w') as f:
        f.write("# PyMOL visualization script generated by PICA\n")
        f.write(f"load {os.path.basename(pdb_path)}\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write(f"color marine, chain {chain_A_id}\n")
        f.write(f"color orange, chain {chain_B_id}\n\n")
        
        # Color interface residues
        f.write("# Select and show interface residues\n")
        res_A_list = "+".join(sorted(list({str(m['res_A_id']) for m in matches})))
        res_B_list = "+".join(sorted(list({str(m['res_B_id']) for m in matches})))
        
        if res_A_list:
            f.write(f"select interface_A, chain {chain_A_id} and resi {res_A_list}\n")
            f.write("color lightblue, interface_A\n")
            f.write("show sticks, interface_A\n")
        if res_B_list:
            f.write(f"select interface_B, chain {chain_B_id} and resi {res_B_list}\n")
            f.write("color lightorange, interface_B\n")
            f.write("show sticks, interface_B\n")
            
        f.write("\n# Represent specific chemical interactions as dashed lines\n")
        f.write("set dash_width, 2.5\n")
        f.write("set dash_color, yellow\n")
        
        counts = {"Hydrophobic": 0, "Salt Bridge": 0, "Hydrogen Bond": 0, "Pi-Pi Stacking": 0}
        for i, m in enumerate(matches):
            itype = m['type']
            counts[itype] = counts.get(itype, 0) + 1
            idx = counts[itype]
            
            # Select atoms for line visualization
            atom_A_sel = f"(chain {chain_A_id} and resi {m['res_A_id']} and name {m['atom_A_name']})"
            atom_B_sel = f"(chain {chain_B_id} and resi {m['res_B_id']} and name {m['atom_B_name']})"
            
            # Use specific color schemes for each bond type
            color_map = {
                "Hydrophobic": "forest",
                "Salt Bridge": "red",
                "Hydrogen Bond": "cyan",
                "Pi-Pi Stacking": "magenta"
            }
            color = color_map.get(itype, "grey")
            
            dist_name = f"dist_{itype.lower().replace(' ', '_')}_{idx}"
            f.write(f"distance {dist_name}, {atom_A_sel}, {atom_B_sel}\n")
            f.write(f"color {color}, {dist_name}\n")
            
        f.write("\nutil.cbc\n")
        f.write("deselect\n")
        print(f"-> PyMOL script successfully generated: {filepath}")


def main():
    # --- COMMAND-LINE INTERFACE CONFIGURATION ---
    parser = argparse.ArgumentParser(
        description="PICA - Protein Interface Contact Analyzer",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pdb_file", help="Path to the PDB file to analyze.")
    parser.add_argument("-a", "--chainA", default="A", help="ID of Chain A (default: 'A').")
    parser.add_argument("-b", "--chainB", default="B", help="ID of Chain B (default: 'B').")
    parser.add_argument("-t", "--threshold", type=float, default=6.0, 
                        help="Global interface distance threshold in Angstroms (default: 6.0 Å).")
    parser.add_argument("-o", "--output", default="contacts_output.csv", 
                        help="Filename for CSV output (default: 'contacts_output.csv').")
    parser.add_argument("-p", "--pymol", default="visualize_contacts.pml", 
                        help="Filename for the PyMOL script (default: 'visualize_contacts.pml').")
    
    args = parser.parse_args()
    
    # Check if the file exists
    if not os.path.exists(args.pdb_file):
        print(f"Error: The PDB file '{args.pdb_file}' could not be found.")
        sys.exit(1)
        
    print("="*80)
    print("                      PROTEIN INTERFACE CONTACT ANALYZER")
    print("="*80)
    print(f"Structure File : {args.pdb_file}")
    print(f"Analyzing      : Chain {args.chainA} <-> Chain {args.chainB}")
    print(f"Interface Cutoff: {args.threshold} Å")
    print("="*80)
    
    # --- PARSING AND CLEANING ---
    try:
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("complex", args.pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file: {e}")
        sys.exit(1)
        
    model = structure[0]
    
    if args.chainA not in model or args.chainB not in model:
        available_chains = [c.get_id() for c in model]
        print(f"Error: Specified chains must be present in the model.")
        print(f"Available chains in this structure: {available_chains}")
        sys.exit(1)
        
    chain_A = model[args.chainA]
    chain_B = model[args.chainB]
    
    atoms_A = get_clean_atoms(chain_A)
    atoms_B = get_clean_atoms(chain_B)
    
    print(f"-> Chain {args.chainA} clean heavy atoms: {len(atoms_A)}")
    print(f"-> Chain {args.chainB} clean heavy atoms: {len(atoms_B)}")
    
    # --- SPATIAL DISTANCE COMPUTATION ---
    print("\n[1] Detecting general spatial contacts...")
    spatial_contacts = []
    for atom_A in atoms_A:
        coord_A = atom_A.get_coord()
        for atom_B in atoms_B:
            coord_B = atom_B.get_coord()
            dist = calculate_distance(coord_A, coord_B)
            if dist <= args.threshold:
                spatial_contacts.append((atom_A, atom_B, dist))
                
    print(f"   Detected {len(spatial_contacts)} raw heavy atom contacts under {args.threshold} Å.")
    
    # --- RECLASSIFYING SPECIFIC INTERESTING NON-COVALENT FORCES ---
    print("\n[2] Filtering and characterizing biochemical interactions...")
    
    detected_matches = []
    
    # 1. Hydrophobic Interactions
    hydrophobic_contacts = []
    for atom_A, atom_B, dist in spatial_contacts:
        res_A = atom_A.get_parent()
        res_B = atom_B.get_parent()
        if res_A.get_resname() in HYDROPHOBIC_RES and res_B.get_resname() in HYDROPHOBIC_RES:
            if atom_A.get_name() not in BACKBONE_ATOMS and atom_B.get_name() not in BACKBONE_ATOMS:
                if atom_A.get_name().startswith('C') and atom_B.get_name().startswith('C'):
                    if dist <= 4.0:
                        hydrophobic_contacts.append((atom_A, atom_B, dist))
                        detected_matches.append({
                            'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                            'res_B_name': res_B.get_resname(), 'res_B_id': res_B.get_id()[1], 'atom_B_name': atom_B.get_name(),
                            'type': 'Hydrophobic', 'distance': dist, 'info': 'Side-chain Carbon-Carbon'
                        })

    # 2. Salt Bridges
    salt_bridges = []
    for atom_A, atom_B, dist in spatial_contacts:
        res_A = atom_A.get_parent()
        res_B = atom_B.get_parent()
        # Direct ionics
        for pos_res, neg_res, cat_atoms, ani_atoms in [
            (res_A, res_B, N_CATION, O_ANION),
            (res_B, res_A, N_CATION, O_ANION)
        ]:
            if pos_res.get_resname() in RES_POSITIFS and neg_res.get_resname() in RES_NEGATIFS:
                # Find matching heavy ions
                if atom_A.get_name() in (cat_atoms if pos_res == res_A else ani_atoms) and \
                   atom_B.get_name() in (ani_atoms if pos_res == res_A else cat_atoms):
                    if dist <= 4.0:
                        salt_bridges.append((atom_A, atom_B, dist))
                        detected_matches.append({
                            'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                            'res_B_name': res_B.get_resname(), 'res_B_id': res_B.get_id()[1], 'atom_B_name': atom_B.get_name(),
                            'type': 'Salt Bridge', 'distance': dist, 'info': 'Charged Cation-Anion contact'
                        })
                        
    # 3. Hydrogen Bonds
    h_bonds = []
    for atom_A, atom_B, dist in spatial_contacts:
        if atom_A.element.strip() in {"N", "O"} and atom_B.element.strip() in {"N", "O"}:
            if dist <= 3.5:
                res_A = atom_A.get_parent()
                ant_A = [atom for atom in res_A if atom != atom_A]
                if ant_A:
                    coord_ant_A = ant_A[0].get_coord()
                    angle = calculate_angle(atom_A.get_coord(), coord_ant_A, atom_B.get_coord())
                    if angle > 90.0:
                        h_bonds.append((atom_A, atom_B, dist, angle))
                        detected_matches.append({
                            'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                            'res_B_name': atom_B.get_parent().get_resname(), 'res_B_id': atom_B.get_parent().get_id()[1], 'atom_B_name': atom_B.get_name(),
                            'type': 'Hydrogen Bond', 'distance': dist, 'info': f'Angle: {angle:.1f}°'
                        })

    # 4. Pi-Pi Stacking (Aromatics)
    aromatic_contacts = []
    res_aromatic_A = [res for res in chain_A if res.get_resname() in AROMATIC_RES and is_res(res)]
    res_aromatic_B = [res for res in chain_B if res.get_resname() in AROMATIC_RES and is_res(res)]
    
    for res_A in res_aromatic_A:
        centroid_A, normal_A = calculate_centroid_normal(res_A)
        if centroid_A is None:
            continue
        for res_B in res_aromatic_B:
            centroid_B, normal_B = calculate_centroid_normal(res_B)
            if centroid_B is None:
                continue
                
            dist_centroids = np.linalg.norm(centroid_A - centroid_B)
            if dist_centroids <= 6.0:
                dot_product = np.dot(normal_A, normal_B)
                cos_angle = min(1.0, max(-1.0, abs(dot_product)))
                angle = np.degrees(np.arccos(cos_angle))
                
                if angle <= 30.0:
                    stack_type = "Parallel (Face-to-Face)"
                elif 60.0 <= angle <= 90.0:
                    stack_type = "T-shaped (Edge-to-Face)"
                else:
                    stack_type = "Displaced/Intermediate"
                    
                aromatic_contacts.append((res_A, res_B, dist_centroids, angle, stack_type))
                
                # To visualize easily in PyMOL, map to their respective CG atom (ring entry)
                atom_A = res_A["CG"] if "CG" in res_A else list(res_A.get_atoms())[0]
                atom_B = res_B["CG"] if "CG" in res_B else list(res_B.get_atoms())[0]
                
                detected_matches.append({
                    'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                    'res_B_name': res_B.get_resname(), 'res_B_id': res_B.get_id()[1], 'atom_B_name': atom_B.get_name(),
                    'type': 'Pi-Pi Stacking', 'distance': dist_centroids, 'info': f'Angle: {angle:.1f}° | Shape: {stack_type}'
                })

    # Print summary statistics to the console
    print("\n" + "="*80)
    print("                           SUMMARY OF DETECTED CONTACTS")
    print("="*80)
    print(f"Hydrophobic Contacts   : {len(hydrophobic_contacts)}")
    print(f"Salt Bridges           : {len(salt_bridges)}")
    print(f"Hydrogen Bonds         : {len(h_bonds)}")
    print(f"Pi-Pi Stacking Contacts: {len(aromatic_contacts)}")
    print("="*80)
    
    # Detailed text log
    if detected_matches:
        print("\nDetail of specific molecular interactions:")
        for m in detected_matches:
            print(f" * [{m['type']}] {m['res_A_name']}{m['res_A_id']} ({args.chainA}) <-> "
                  f"{m['res_B_name']}{m['res_B_id']} ({args.chainB}) | Dist: {m['distance']:.2f} Å | {m['info']}")
    else:
        print("\nNo specific non-covalent interactions found at the strict cutoffs.")
    print("="*80)
    
    # --- WRITING EXPORTS ---
    # CSV generation
    try:
        with open(args.output, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(["Chain_A_Resname", "Chain_A_ResID", "Chain_A_Atom", 
                             "Chain_B_Resname", "Chain_B_ResID", "Chain_B_Atom", 
                             "Interaction_Type", "Distance_Angstrom", "Details"])
            for m in detected_matches:
                writer.writerow([m['res_A_name'], m['res_A_id'], m['atom_A_name'],
                                 m['res_B_name'], m['res_B_id'], m['atom_B_name'],
                                 m['type'], round(m['distance'], 3), m['info']])
        print(f"-> CSV data report successfully written: {args.output}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")
        
    # PyMOL PML script generation
    write_pymol_script(args.pymol, args.pdb_file, args.chainA, args.chainB, detected_matches)
    print("="*80)
    print("Analysis complete. You can now load the generated .pml script directly in PyMOL!")
    print("="*80)


if __name__ == "__main__":
    main()
