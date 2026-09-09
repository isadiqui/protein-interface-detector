#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Detection and Classification of Protein Interface Contacts.

This script coordinates the parsing of PDB files, executes distance-based search
algorithms, filters and classifies non-covalent interactions( hydrophobic, ionic,
hydrogen bonds, pi-pi stacking) at the interface of 2 chains, exports
results to CSV reports and PyMOL visualization scripts.

"""
import os
import sys
import csv
import argparse
import numpy as np
from Bio.PDB import PDBParser




# we import local modules
from biochem_geometry import ( 
    HYDROPHOBIC_RES, BACKBONE_ATOMS, RES_POSITIFS, RES_NEGATIFS,
    N_CATION, O_ANION, AROMATIC_RES,
    get_clean_atoms, calculate_distance, calculate_angle,
    calculate_centroid_normal, is_res
)

from pymol_generator import write_pymol_script



def main():
    """
    Main execution function for the analysis pipeline.

    Configures command-line arguments, parses structural coordinates, identifies
    all atomic paires matching distance thresholds, filters specific chemical forces
    and handles the production of text console reports, CSV outputs and PyMOL visualization.

    Raises:
        SystemExit: If the PDB file does not exist, structural parsing fails or
        the specified chains are absent from the structure.

    """
    parser = argparse.ArgumentParser(
        description="Protein Interface Contact Finder",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("-a", "--chainA", default="A", help="ID of chain A (default: 'A')")
    parser.add_argument("-b", "--chainB", default="B", help="ID of chain B (default: 'B')")
    parser.add_argument("-t", "--threshold", type=float, default=6.0,
                        help="Global interface distance threshold in Angstroms (default: 6.0 A)")
    parser.add_argument("-o", "--output", default="results/contacts_output.csv",
                        help="Filename for CSV output (default: 'contacts_output.csv')")
    parser.add_argument("-p", "--pymol", default="results/contacts.pml",
                        help="Filename for the Pymol script (default: 'contacts.pml')")
    args = parser.parse_args()

    # Check if file exists 
    if not os.path.exists(args.pdb_file):
        print(f"Error: The PDB file '{args.pdb_file}' could not be found")
        sys.exit(1)

    print(f"Structure file: {args.pdb_file}")
    print(f"Analyzing chain {args.chainA} and chain {args.chainB}")
    print(f"Cutoff: {args.threshold}")
    print("\n")

    # Parsing and cleaning
    try:
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("complex", args.pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file: {e}")
        sys.exit(1)

    model = structure[0]

    if args.chainA not in model or args.chainB not in model:
        available_chains = [c.get_id() for c in model]
        print(f"Error: Specified chains must be present in the model")
        print(f"Available chains in this structure: {available_chains}")
        sys.exit(1)

    chain_A = model[args.chainA]
    chain_B = model[args.chainB]

    atoms_A = get_clean_atoms(chain_A)
    atoms_B = get_clean_atoms(chain_B)

    print(f"Chain {args.chainA} clean heavy atoms: {len(atoms_A)}")
    print(f"Chain {args.chainB} clean heavy atoms: {len(atoms_B)}")
    print("\n")

    # Detecting general contacts
    contacts = []
    for atom_A in atoms_A:
        coord_A = atom_A.get_coord()
        for atom_B in atoms_B:
            coord_B = atom_B.get_coord()
            dist = calculate_distance(coord_A, coord_B)
            if dist <= args.threshold:
                contacts.append((atom_A, atom_B, dist))
    print(f"Detected {len(contacts)} raw heavy atom contacts under {args.threshold} Å")
    print("\n")

    # Filtering and characterizing biochemical interactions
    matches = []

    # Hydrophobic interactions
    hydrophobic_contacts = []
    for atom_A, atom_B, dist in contacts:
        res_A = atom_A.get_parent()
        res_B = atom_B.get_parent()
        if res_A.get_resname() in HYDROPHOBIC_RES and res_B.get_resname() in HYDROPHOBIC_RES:
            if atom_A.get_name() not in BACKBONE_ATOMS and atom_B.get_name() not in BACKBONE_ATOMS:
                if atom_A.get_name().startswith('C') and atom_B.get_name().startswith('C'):
                    if dist <= 4.0:
                        hydrophobic_contacts.append((atom_A, atom_B, dist))
                        matches.append({
                            'res_A_name' : res_A.get_resname(), 'res_A_id' : res_A.get_id()[1], 'atom_A_name' : atom_A.get_name(),
                            'res_B_name' : res_B.get_resname(), 'res_B_id' : res_B.get_id()[1], 'atom_B_name' : atom_B.get_name(),
                            'type': 'Hydrophobic', 'distance': dist})

    # Salt bridges
    salt_bridges = []
    for atom_A, atom_B, dist in contacts:
        res_A = atom_A.get_parent()
        res_B = atom_B.get_parent()
        # direct ionics
        for pos_res, neg_res, cat_atoms, ani_atoms in [
            (res_A, res_B, N_CATION, O_ANION),
            (res_B, res_A, N_CATION, O_ANION)
        ]:
            if pos_res.get_resname() in RES_POSITIFS and neg_res.get_resname() in RES_NEGATIFS:
                if atom_A.get_name() in (cat_atoms if pos_res == res_A else ani_atoms) and \
                   atom_B.get_name() in (ani_atoms if pos_res == res_A else cat_atoms):
                   if dist <= 4.0:
                       salt_bridges.append((atom_A, atom_B, dist))
                       matches.append({
                           'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                           'res_B_name': res_B.get_resname(), 'res_B_id': res_B.get_id()[1], 'atom_B_name': atom_B.get_name(),
                           'type': 'Salt bridge', 'distance': dist
                       })

    # Hydrogen bonds
    h_bonds = []
    for atom_A, atom_B, dist in contacts:
        if atom_A.element.strip() in {"N", "O"} and atom_B.element.strip() in {"N", "O"}:
            if dist <= 3.5:
                res_A = atom_A.get_parent()
                ant_A = [atom for atom in res_A if atom != atom_A]
                if ant_A:
                    coord_ant_A = ant_A[0].get_coord()
                    angle = calculate_angle(atom_A.get_coord(),coord_ant_A, atom_B.get_coord())
                    if angle > 90.0:
                        h_bonds.append((atom_A, atom_B, dist, angle))
                        matches.append({
                            'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                            'res_B_name': atom_B.get_parent().get_resname(), 'res_B_id': atom_B.get_parent().get_id()[1], 'atom_B_name': atom_B.get_name(),
                            'type': 'Hydrogen bond', 'distance': dist
                        })

    # Pi-Pi stacking
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
                    stack_type = "Parallel"
                elif 60.0 <= angle <= 90.0:
                    stack_type = "T-shaped"
                else:
                    stack_type = "Displaced/Intermediate"

                aromatic_contacts.append((res_A, res_B, dist_centroids, angle, stack_type))

                # To visualize in PyMOL
                atom_A = res_A["CG"] if "CG" in res_A else list(res_A.get_atoms())[0]
                atom_B = res_B["CG"] if "CG" in res_B else list(res_B.get_atoms())[0]

                matches.append({
                    'res_A_name': res_A.get_resname(), 'res_A_id': res_A.get_id()[1], 'atom_A_name': atom_A.get_name(),
                    'res_B_name': res_B.get_resname(), 'res_B_id': res_B.get_id()[1], 'atom_B_name': atom_B.get_name(),    
                    'type': 'Pi-Pi stacking', 'distance': dist_centroids
                })

    # Summary 
    print("      SUMMARY OF DETECTED CONTACTS")
    print(f"Hydrophobic contacts   : {len(hydrophobic_contacts)}")
    print(f"Salt bridges           : {len(salt_bridges)}")
    print(f"Hydrogen bonds         : {len(h_bonds)}")
    print(f"Pi-Pi stacking contacts: {len(aromatic_contacts)}")
    print("\n")

    if matches:
        print("\n DETAILS:")
        for m in matches:
            print(f" [{m['type']}] {m['res_A_name']}{m['res_A_id']} ({args.chainA}) <->"
                  f"{m['res_B_name']}{m['res_B_id']} ({args.chainB}) | Distance: {m['distance']:.2f} Å")
    print("\n")  

    # Exports
    # Generate CSV
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    pymol_dir = os.path.dirname(args.pymol)
    if pymol_dir:
        os.makedirs(pymol_dir, exist_ok=True)

    try: 
        with open(args.output, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(["Chain_A_Resname", "Chain_A_ResID", "Chain_A_Atom",
                             "Chain_B_Resname", "Chain_B_ResID", "Chain_B_Atom",
                             "Interaction_Type", "Distance_Angstrom", "Details"])
            for m in matches:
                writer.writerow([m['res_A_name'], m["res_A_id"], m['atom_A_name'],
                                 m['res_B_name'], m["res_B_id"], m['atom_B_name'],
                                 m['type'], round(m['distance'], 3)])
            print(f"CSV data report: {args.output}")

    except Exception as e:
        print(f"Error writing CSV file {e}")

    #PyMOL script
    write_pymol_script(args.pymol, args.pdb_file, args.chainA, args.chainB, matches)

if __name__ == "__main__":
    main()