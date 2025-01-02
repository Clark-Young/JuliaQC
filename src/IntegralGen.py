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
import numpy as np
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
    basis_set = rest[0]
    
    atoms_lines = lines[2:2 + num_atoms]
    atoms = ''.join(atoms_lines).strip()
    
    return charge, multiplicity, basis_set, atoms


if len(sys.argv) != 2:
    print("Usage: python IntegralGen.py input.xyz")
    sys.exit(1)

input_file = sys.argv[1]
charge, multiplicity,basis_set, atoms = read_xyz_file(input_file)

mol = pyscf.gto.M(
    atom = atoms,
    basis = basis_set,
    charge = charge,
    spin = multiplicity - 1,
    verbose = 0,
)

mf = scf.RHF(mol)
rhf_obj = pyscf.scf.RHF(mol)
rhf_obj.max_cycle = 200
rhf_obj.conv_tol  = 1e-10
dm_rhf = rhf_obj.get_init_guess()
rhf_obj.kernel(dm_rhf) 
if rhf_obj.converged:
    hcore  = mol.intor('int1e_nuc')
    hcore += mol.intor('int1e_kin')
    ovlp   = mol.intor('int1e_ovlp')
    eri    = mol.intor('int2e')

    # 保存积分
    np.save('h_core.npy', hcore)
    np.save('eri.npy', eri)
    np.save('ovlp.npy', ovlp)

else:
    print("RHF calculation failed.")
    exit(1)





