# Description: This Python Code is used to generate the integrals for the quantum chemistry calculation later.
#   The integrals include overlap, kinetic, nuclear attraction, electron repulsion, and two electron integrals.
#   It's very troublesome to do the double electron integral calculation, so I use the PySCF package to do this currently.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-13
# Reference: Junjie Yang's code on GitHub, https://github.com/yangjunjie0320/hf-tutorial


import os
import sys
import numpy
import pyscf
from pyscf import gto
from pyscf import scf, fci

def read_xyz_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    num_atoms = int(lines[0].strip())
    
    charge, multiplicity, *rest = lines[1].split()
    charge = int(charge)
    multiplicity = int(multiplicity)
    
    atoms_lines = lines[2:2 + num_atoms]
    atoms = ''.join(atoms_lines).strip()
    
    return charge, multiplicity, atoms


if len(sys.argv) != 2:
    print("Usage: python IntegralGen.py input.xyz")
    sys.exit(1)

input_file = sys.argv[1]
charge, multiplicity, atoms = read_xyz_file(input_file)





