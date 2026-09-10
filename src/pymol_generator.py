#-*- coding: utf-8 -*-
"""
PyMOL script generator module.

This module provides automated generation of PyMOL visualization files,
enabling immediate 3D representation of molecular interactions at the interface.

"""
import os

def write_pymol_script(filepath, pdb_path, chain_A_id, chain_B_id, matches):
    """ 
    Generates a customized PyMOL script to visualize detected interactions.

    Writes command lines to load the PDB structure, color the respective chains, 
    display interface residues as sticks and represent each chemical interaction 
    using specific dashed line coloring schemes.

    Args:
        filepath (str): Path where the output PyMOL script will be saved.
        pdb_path (str): Path to the source PDB file being analyzed.
        chain_A_id (str): Chain ID of the receptor protein.
        chain_B_id (str): Chain ID of the binding partner.
        matches (list of dict): List of detected interaction dictionaries. Each 
        dictionary must contain keys like 'res_A_id', 'res_B_id', 'atom_A_name', 
        'atom_B_name', 'type', and 'distance'.

    Returns : 
        None

    """

    with open(filepath, 'w') as f:
        abs_pdb_path = os.path.abspath(pdb_path).replace("\\", "/")
        if abs_pdb_path.startswith("/mnt/"):
            parts = abs_pdb_path.split("/")
            drive_letter = parts[2].upper()
            abs_pdb_path = f"{drive_letter}:/" + "/".join(parts[3:])

        f.write(f"load {abs_pdb_path}\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write(f"color blue, chain {chain_A_id}\n")
        f.write(f"color red, chain {chain_B_id}\n\n")

        # Color interface residues
        f.write("#select and show interface residues\n")
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

    
        counts = {"Hydrophobic": 0, "Salt bridge": 0, "Hydrogen bond": 0, "Pi-Pi stacking": 0}
        for i, m in enumerate(matches):
            itype = m['type']
            counts[itype] = counts.get(itype, 0) + 1
            idx = counts[itype]

            # Select atoms for line visualisation
            atom_A_sel = f"(chain {chain_A_id} and resi {m['res_A_id']} and name {m['atom_A_name']})"
            atom_B_sel = f"(chain {chain_B_id} and resi {m['res_B_id']} and name {m['atom_B_name']})"

            # Different colors for different bond types
            color_map = {
                "Hydrophobic": "forest",
                "Salt bridge": "red", 
                "Hydrogen bond": "cyan", 
                "Pi-Pi stacking": "magenta"
            }
            color = color_map.get(itype, "grey")

            dist_name = f"dist_{itype.lower().replace(' ', '_')}_{idx}"
            f.write(f"distance {dist_name}, {atom_A_sel}, {atom_B_sel}\n")
            f.write(f"color {color}, {dist_name}\n")
            f.write(f"set dash_color, {color}, {dist_name}\n")

        f.write("\nutil.cnc\n")
        f.write("deselect\n")
        print(f"PyMOL script generated: {filepath}")