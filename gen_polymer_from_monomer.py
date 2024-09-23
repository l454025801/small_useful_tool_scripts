#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:46:19 2022

@author: leon

Create an inp file for packmol to generate a pdb file of a long polymer chain. The chain will extend along x axis. PDB file for single unit is required.
The relative coordinates of 3 atoms in one repeat unit are required to fix the configuration. The length of the repeat unit is need to extend the chain.

"""

def write_inp(chain_length, atom_index, atom_coordinates, unit_length, repeat_unit_name):
    '''
    Parameters
    ----------
    chain_length : int
                    Give a integer and generate an .inp file to write a polymer pdb file.
    atom_index   : list of int
                    indices of 3 atoms used to fix the coordinates. Please use 2 connecting beckbone atoms for atom 1 and atom 2
    atom_coordinates : list of lists 
                    contain the coordinates of 3 atoms of atom_index; eg: [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]
                    recommand to have x1=y1=z1=0 and all z=0
    unit_length  : float
                    the length of the repeat unit, PLUS the length of the bond connecting with the next repeat unit
    repeat_unit_name : str
                    the file name of the repeat unit pdb file

    Returns
    -------
    None.

    '''

    x1,y1,z1 = atom_coordinates[0]
    x2,y2,z2 = atom_coordinates[1]
    x3,y3,z3 = atom_coordinates[2]

    if type(chain_length) != int:
        return ('Wrong chain length')
    
    filename = "build_long_chain_" + str(chain_length) +".inp"
    inp = open(filename, "w") 
    inp_script = \
    "tolerance 0.5\n\
    output polymer_DOP_%i.pdb\n\
    filetype pdb\n\
    seed -1\n\
    \n\
        " %(chain_length)
    
    for i in range(chain_length):
        # 
        # if the number of beckbone atoms is even, then you do not need this if statement, 
        # otherwise need to consider that the starting atom of second unit is not at the same y,z position as it of the first unit
        # eg:
        #   1   3   2          1   3   1   3
        #    \ / \ / \    or    \ / \ / \ / \
        #     2   1   3          2   4   2   4
        # This function provides a sloopy way to deal this condition by simply vertically flipping the chain. 
        # You probably need to equilibrate the chain before actually use it.
        if i%2 == 0:
            inp_script += \
                    "structure %s\n\
     number 1 \n\
     chain A \n\
     resnumbers 2 \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     end structure \n\
        \n\
         " %(repeat_unit_name, atom_index[0], x1+i*unit_length, y1, z1, atom_index[1], x2+i*unit_length, y2, z2, atom_index[2], x3+i*unit_length, y3, z3)
        else:
            inp_script += \
                    "structure %s\n\
     number 1 \n\
     chain A\n\
     resnumbers 2 \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     atoms %i \n\
     inside sphere %.2f %.2f %.2f 0.05 \n\
     end atoms \n\
     end structure \n\
         \n\
         " %(repeat_unit_name, atom_index[0], x1+i*unit_length, -y1, -z1, atom_index[1], x2+i*unit_length, -y2, -z2, atom_index[2], x3+i*unit_length, -y3, -z3)

    inp.write(inp_script)    
        
    print ("finish writing .inp file")
    return None

