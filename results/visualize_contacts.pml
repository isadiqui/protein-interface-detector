load C:/Users/hp/projet1/data/2XA0.pdb
hide everything
show cartoon
color blue, chain A
color red, chain C

#select and show interface residues
select interface_A, chain A and resi 104+107+108+110+112+118+133+136+137+139+140+143+144+145+146+153+200+201+204
color lightblue, interface_A
show sticks, interface_A
select interface_B, chain C and resi 57+59+61+63+64+66+68+69+70+71+73+74+78+81
color lightorange, interface_B
show sticks, interface_B
distance dist_hydrophobic_1, (chain A and resi 104 and name CB), (chain C and resi 70 and name CD1)
color forest, dist_hydrophobic_1
set dash_color, forest, dist_hydrophobic_1
distance dist_hydrophobic_2, (chain A and resi 104 and name CG), (chain C and resi 66 and name CG2)
color forest, dist_hydrophobic_2
set dash_color, forest, dist_hydrophobic_2
distance dist_hydrophobic_3, (chain A and resi 104 and name CD1), (chain C and resi 66 and name CG2)
color forest, dist_hydrophobic_3
set dash_color, forest, dist_hydrophobic_3
distance dist_hydrophobic_4, (chain A and resi 104 and name CE1), (chain C and resi 66 and name CB)
color forest, dist_hydrophobic_4
set dash_color, forest, dist_hydrophobic_4
distance dist_hydrophobic_5, (chain A and resi 104 and name CE1), (chain C and resi 66 and name CG2)
color forest, dist_hydrophobic_5
set dash_color, forest, dist_hydrophobic_5
distance dist_hydrophobic_6, (chain A and resi 104 and name CZ), (chain C and resi 66 and name CB)
color forest, dist_hydrophobic_6
set dash_color, forest, dist_hydrophobic_6
distance dist_hydrophobic_7, (chain A and resi 104 and name CZ), (chain C and resi 66 and name CG2)
color forest, dist_hydrophobic_7
set dash_color, forest, dist_hydrophobic_7
distance dist_hydrophobic_8, (chain A and resi 108 and name CD2), (chain C and resi 66 and name CG1)
color forest, dist_hydrophobic_8
set dash_color, forest, dist_hydrophobic_8
distance dist_hydrophobic_9, (chain A and resi 108 and name CE2), (chain C and resi 63 and name CD1)
color forest, dist_hydrophobic_9
set dash_color, forest, dist_hydrophobic_9
distance dist_hydrophobic_10, (chain A and resi 112 and name CE1), (chain C and resi 66 and name CG1)
color forest, dist_hydrophobic_10
set dash_color, forest, dist_hydrophobic_10
distance dist_hydrophobic_11, (chain A and resi 112 and name CE1), (chain C and resi 66 and name CD1)
color forest, dist_hydrophobic_11
set dash_color, forest, dist_hydrophobic_11
distance dist_hydrophobic_12, (chain A and resi 133 and name CG1), (chain C and resi 59 and name CD2)
color forest, dist_hydrophobic_12
set dash_color, forest, dist_hydrophobic_12
distance dist_hydrophobic_13, (chain A and resi 133 and name CG1), (chain C and resi 63 and name CD2)
color forest, dist_hydrophobic_13
set dash_color, forest, dist_hydrophobic_13
distance dist_hydrophobic_14, (chain A and resi 144 and name CB), (chain C and resi 74 and name CE)
color forest, dist_hydrophobic_14
set dash_color, forest, dist_hydrophobic_14
distance dist_hydrophobic_15, (chain A and resi 153 and name CB), (chain C and resi 63 and name CD2)
color forest, dist_hydrophobic_15
set dash_color, forest, dist_hydrophobic_15
distance dist_hydrophobic_16, (chain A and resi 153 and name CG), (chain C and resi 63 and name CD2)
color forest, dist_hydrophobic_16
set dash_color, forest, dist_hydrophobic_16
distance dist_hydrophobic_17, (chain A and resi 153 and name CE1), (chain C and resi 59 and name CD2)
color forest, dist_hydrophobic_17
set dash_color, forest, dist_hydrophobic_17
distance dist_hydrophobic_18, (chain A and resi 201 and name CD2), (chain C and resi 74 and name CB)
color forest, dist_hydrophobic_18
set dash_color, forest, dist_hydrophobic_18
distance dist_hydrophobic_19, (chain A and resi 201 and name CD2), (chain C and resi 74 and name CG)
color forest, dist_hydrophobic_19
set dash_color, forest, dist_hydrophobic_19
distance dist_hydrophobic_20, (chain A and resi 204 and name CB), (chain C and resi 81 and name CB)
color forest, dist_hydrophobic_20
set dash_color, forest, dist_hydrophobic_20
distance dist_salt_bridge_1, (chain A and resi 107 and name NH1), (chain C and resi 69 and name OE1)
color red, dist_salt_bridge_1
set dash_color, red, dist_salt_bridge_1
distance dist_salt_bridge_2, (chain A and resi 110 and name NH1), (chain C and resi 69 and name OE1)
color red, dist_salt_bridge_2
set dash_color, red, dist_salt_bridge_2
distance dist_salt_bridge_3, (chain A and resi 110 and name NH1), (chain C and resi 69 and name OE2)
color red, dist_salt_bridge_3
set dash_color, red, dist_salt_bridge_3
distance dist_salt_bridge_4, (chain A and resi 110 and name NH2), (chain C and resi 69 and name OE1)
color red, dist_salt_bridge_4
set dash_color, red, dist_salt_bridge_4
distance dist_salt_bridge_5, (chain A and resi 110 and name NH2), (chain C and resi 69 and name OE2)
color red, dist_salt_bridge_5
set dash_color, red, dist_salt_bridge_5
distance dist_salt_bridge_6, (chain A and resi 139 and name NH2), (chain C and resi 61 and name OE1)
color red, dist_salt_bridge_6
set dash_color, red, dist_salt_bridge_6
distance dist_salt_bridge_7, (chain A and resi 140 and name OD1), (chain C and resi 64 and name NH2)
color red, dist_salt_bridge_7
set dash_color, red, dist_salt_bridge_7
distance dist_salt_bridge_8, (chain A and resi 140 and name OD2), (chain C and resi 64 and name NH2)
color red, dist_salt_bridge_8
set dash_color, red, dist_salt_bridge_8
distance dist_salt_bridge_9, (chain A and resi 146 and name NH1), (chain C and resi 68 and name OD1)
color red, dist_salt_bridge_9
set dash_color, red, dist_salt_bridge_9
distance dist_salt_bridge_10, (chain A and resi 146 and name NH1), (chain C and resi 68 and name OD2)
color red, dist_salt_bridge_10
set dash_color, red, dist_salt_bridge_10
distance dist_salt_bridge_11, (chain A and resi 200 and name OE2), (chain C and resi 78 and name NH2)
color red, dist_salt_bridge_11
set dash_color, red, dist_salt_bridge_11
distance dist_hydrogen_bond_1, (chain A and resi 107 and name NH2), (chain C and resi 73 and name ND2)
color cyan, dist_hydrogen_bond_1
set dash_color, cyan, dist_hydrogen_bond_1
distance dist_hydrogen_bond_2, (chain A and resi 110 and name NH1), (chain C and resi 69 and name OE2)
color cyan, dist_hydrogen_bond_2
set dash_color, cyan, dist_hydrogen_bond_2
distance dist_hydrogen_bond_3, (chain A and resi 110 and name NH2), (chain C and resi 69 and name OE1)
color cyan, dist_hydrogen_bond_3
set dash_color, cyan, dist_hydrogen_bond_3
distance dist_hydrogen_bond_4, (chain A and resi 118 and name OE1), (chain C and resi 57 and name O)
color cyan, dist_hydrogen_bond_4
set dash_color, cyan, dist_hydrogen_bond_4
distance dist_hydrogen_bond_5, (chain A and resi 118 and name OE1), (chain C and resi 59 and name N)
color cyan, dist_hydrogen_bond_5
set dash_color, cyan, dist_hydrogen_bond_5
distance dist_hydrogen_bond_6, (chain A and resi 136 and name O), (chain C and resi 64 and name NE)
color cyan, dist_hydrogen_bond_6
set dash_color, cyan, dist_hydrogen_bond_6
distance dist_hydrogen_bond_7, (chain A and resi 137 and name O), (chain C and resi 64 and name NE)
color cyan, dist_hydrogen_bond_7
set dash_color, cyan, dist_hydrogen_bond_7
distance dist_hydrogen_bond_8, (chain A and resi 139 and name NE), (chain C and resi 64 and name NH1)
color cyan, dist_hydrogen_bond_8
set dash_color, cyan, dist_hydrogen_bond_8
distance dist_hydrogen_bond_9, (chain A and resi 140 and name OD1), (chain C and resi 64 and name NH2)
color cyan, dist_hydrogen_bond_9
set dash_color, cyan, dist_hydrogen_bond_9
distance dist_hydrogen_bond_10, (chain A and resi 143 and name OD1), (chain C and resi 71 and name OD2)
color cyan, dist_hydrogen_bond_10
set dash_color, cyan, dist_hydrogen_bond_10
distance dist_hydrogen_bond_11, (chain A and resi 143 and name ND2), (chain C and resi 68 and name OD1)
color cyan, dist_hydrogen_bond_11
set dash_color, cyan, dist_hydrogen_bond_11
distance dist_hydrogen_bond_12, (chain A and resi 144 and name N), (chain C and resi 71 and name OD2)
color cyan, dist_hydrogen_bond_12
set dash_color, cyan, dist_hydrogen_bond_12
distance dist_hydrogen_bond_13, (chain A and resi 145 and name N), (chain C and resi 71 and name OD2)
color cyan, dist_hydrogen_bond_13
set dash_color, cyan, dist_hydrogen_bond_13
distance dist_hydrogen_bond_14, (chain A and resi 146 and name NH1), (chain C and resi 68 and name OD2)
color cyan, dist_hydrogen_bond_14
set dash_color, cyan, dist_hydrogen_bond_14
distance dist_hydrogen_bond_15, (chain A and resi 146 and name NH2), (chain C and resi 64 and name NE)
color cyan, dist_hydrogen_bond_15
set dash_color, cyan, dist_hydrogen_bond_15
distance dist_hydrogen_bond_16, (chain A and resi 146 and name NH2), (chain C and resi 64 and name NH2)
color cyan, dist_hydrogen_bond_16
set dash_color, cyan, dist_hydrogen_bond_16
distance dist_hydrogen_bond_17, (chain A and resi 200 and name OE2), (chain C and resi 78 and name NH2)
color cyan, dist_hydrogen_bond_17
set dash_color, cyan, dist_hydrogen_bond_17

util.cnc
deselect
