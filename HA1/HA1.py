'''
Albin Annér, analbin@student.chalmers.se
David Kastö, kasto@student.chalmers.se
'''

#!/usr/bin/env python

# # imports # #
import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import fmin
from ase import Atoms
from ase.calculators.lj import LennardJones

# # functions # #
def V_LJ(r12):
    ''' Function that calculates the Lennard-Jones (LJ) potential between two atoms'''
    global epsilon
    global sigma
    return 4*epsilon*((sigma/r12)**12 - (sigma/r12)**6)

def tot_energy(pos_of_atoms):
    ''' Function that calculates the total energy in a system
        based on an interatomic Lennard-Jones potential '''
    (ni, _) = np.shape(pos_of_atoms) # ni is nbr. of atoms in data file
    total_energy = 0
    # Loops through pair-wise interactions. Interaction between each pair of
    # atoms is only contributing once
    for i in range(ni):
        for j in range(i+1, ni):
            atom_1_pos = pos_of_atoms[i, :]
            atom_2_pos = pos_of_atoms[j,:]
            r12 = np.sqrt(np.sum((atom_1_pos-atom_2_pos)**2))
            total_energy += V_LJ(r12)
    return total_energy

# # input file # #
filename1 = 'coord.txt'  # file with cartesian coordinates for the atoms

# # import and manage data # #
pos_of_atoms = np.array(np.loadtxt(filename1)) #
(ni, nj) = np.shape(pos_of_atoms)

# # Constants # #
epsilon = 0.0104  # eV, parameter for Lennard-Jones potential for Argon
sigma = 3.40  # Å, parameter for Lennard-Jones potential for Argon

# # Problem 1-3  # #
atom_1_pos = np.copy(pos_of_atoms[0,:]) # initial values for optimization procedure
atom_2_pos = np.copy(pos_of_atoms[1,:]) # initial values for optimization procedure
r12 = np.sqrt(np.sum((atom_1_pos-atom_2_pos)**2))
R_0 = fmin(V_LJ, r12) # equilibrium distance
E_own = tot_energy(pos_of_atoms)
print('Equilibrium distance between Ar atoms, assuming LJ potential: ' + str(R_0[0]))
print('Total energy, calculated with our own LJ function: ' + str(E_own))

# # Problem 4 # #
atoms = Atoms('Ar' + str(ni), positions=pos_of_atoms)
calc = LennardJones(epsilon=epsilon,sigma=sigma, rc=100) #create calculator
atoms.set_calculator(calc) #attach calc to Atoms object
E_ASE = atoms.get_potential_energy() #calculate energy
print('Total energy, calculated with ASE: ' + str(E_ASE))

perc_diff=((E_ASE-E_own)/E_ASE)*100
print('Percental difference between our method and the ASE method: ' + str(perc_diff))
