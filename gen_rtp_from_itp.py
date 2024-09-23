#!/usr/bin/env python

# used to create rtp entries for charmm36 forcefield with CGenFF generated itp files.

import numpy as np
import pandas as pd
import os
import argparse


def read_itp(itp_loc):
    print ('Reading itp file')
    itp = open(itp_loc,'r')
    atoms = dict()
    bonds = dict()
    impro = dict()
    copy_atoms = False
    copy_bonds = False
    copy_impro = False
    bond_num = 1
    impro_num = 1
    for line in itp:
        if line.strip() == '[ atoms ]':
            copy_atoms = True
            continue
        elif line.strip() == '[ bonds ]':
            copy_atoms = False
            copy_bonds = True
            continue
        elif line.strip() == '[ pairs ]':
            copy_bonds = False
            continue
        elif line.strip() == '[ dihedrals ]':
            copy_impro = True
            continue
        elif copy_atoms == True:
            if line != '\n' and line.split()[-1] == ';' :
                atoms[str(line.split()[0])] = line.split()
        elif copy_bonds == True:
            if line != '\n' and line.split()[-1] == '1' :
                bonds[str(bond_num)] = line.split()
                bond_num += 1
        elif copy_impro == True:
            if line != '\n' and line.split()[-1] == '2':
                impro[str(impro_num)] = line.split()
                impro_num += 1
    itp.close()
    atoms = pd.DataFrame.from_dict(atoms, orient='index', columns=['#','type','resnum','resname','name','seg_num','charge','mass',';'])
    bonds = pd.DataFrame.from_dict(bonds, orient='index', columns=['atom1','atom2','type'])
    impro = pd.DataFrame.from_dict(impro, orient='index', columns=['atom1','atom2','atom3','atom4','type'])
    print ('------- Done reading itp file -------')
    return atoms, bonds, impro




def write_rtp_entry(atoms,bonds,impro):
    print ('Writing rtp entry')
    resname = atoms.iloc[0]['resname']
    rtp_file = open('rtp_entry.txt','w')
    rtp_file.write('[ %s ]\n' % resname)
    rtp_file.write('; any comments you want \n')
    rtp_file.write('  [ atoms ]\n')
    for i in range(len(atoms)):
        rtp_file.write('      %s   %s   %.4f   %i \n' %(atoms.iloc[i]['name'],atoms.iloc[i]['type'],float(atoms.iloc[i]['charge']),int(atoms.iloc[i]['seg_num'])))
    rtp_file.write('  [ bonds ]\n')
    for i in range(len(bonds)):
        atom1_name = atoms.iloc[int(bonds.iloc[i]['atom1'])-1]['name']
        atom2_name = atoms.iloc[int(bonds.iloc[i]['atom2'])-1]['name']
        rtp_file.write('        %s   %s \n' %(atom1_name,atom2_name))
    rtp_file.write('  [ impropers ]\n')
    for i in range(len(impro)):
        atom1_name = atoms.iloc[int(impro.iloc[i]['atom1'])-1]['name']
        atom2_name = atoms.iloc[int(impro.iloc[i]['atom2'])-1]['name']
        atom3_name = atoms.iloc[int(impro.iloc[i]['atom3'])-1]['name']
        atom4_name = atoms.iloc[int(impro.iloc[i]['atom4'])-1]['name']
        rtp_file.write('        %s   %s    %s   %s \n' %(atom1_name,atom2_name,atom3_name,atom4_name))
    rtp_file.close()
    print ('---------- Done ----------')
    return None




if __name__ == "__main__":
    print ('This is a simple script help generating rtp entries for gromacs\n')
    print ('The connectivity is not considered in this script (ex. +C, -N)\n')
    print ('CMAP is not considered')
    print ('--------------------------------------------------------------------')
    parser = argparse.ArgumentParser()
    parser.add_argument('--itp', type=str, default=False, help='the location of the itp file')
    args = parser.parse_args()
    itp_loc = args.itp
    
    atoms,bonds,impro = read_itp(itp_loc)
    write_rtp_entry(atoms,bonds,impro)






