# Detection and Classification of Protein Interface Contacts

A Python command-line tool to analyze, characterize, and visualize non-covalent interactions at protein-protein interfaces.

## Overview

This project provides an automated pipeline to clean PDB structural data and classify four major types of inter-chain non-covalent interactions:
* **Hydrophobic Contacts:** Carbon-carbon side-chain contacts $\le 4.0\$ Å.
* **Salt Bridges:** Electrostatic contacts $\le 4.0\text{ \AA}$ between cationic nitrogens (Arg, Lys, His) and acidic oxygen anions (Asp, Glu).
* **Hydrogen Bonds:** Nitrogen-oxygen contacts $\le 3.5\text{ \AA}$ validated by an antecedent angle criterion ($\theta > 90.0^\circ$).
* **Pi-Pi Stacking:** Aromatic ring center contacts $\le 6.0\text{ \AA}$ (Phe, Tyr, Trp) categorized into parallel, T-shaped, or displaced geometry.

## Install

Install the required dependencies:

```bash
pip install biopython numpy
```

## Usage

Command-Line Interface (CLI)

Run the script from your terminal:

```bash
python src/detect_contacts.py <pdb_file_path> -a <chain_A> -b <chain_B> [options]
```

Options:

    pdb_file : (Required) Path to the input PDB file (e.g., data/2XA0.pdb).
    -a, --chainA : Chain ID of the first protein partner (default: A).
    -b, --chainB : Chain ID of the second protein partner (default: B).
    -t, --threshold : Maximum interface cutoff distance in Angstroms (default: 6.0 Å).
    -o, --output : Destination path for the CSV report (default: results/contacts_output.csv).
    -p, --pymol : Destination path for the PyMOL script (default: results/visualize_contacts.pml).

Example:

```bash
python src/detect_contacts.py data/2XA0.pdb -a A -b C -o results/contacts.csv -p results/visualize_contacts.pml
```

## PyMOL Visualization

Launch PyMOL.

Open File -> Run Script... and select results/visualize_contacts.pml (or run @results/visualize_contacts.pml in the PyMOL prompt).

Interface residues display as sticks and interactions appear as dashed lines.
   

