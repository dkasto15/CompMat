#!/usr/bin/env python
# coding=utf-8

# # # imports # # #
import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.io import read
from ase.units import kJ, J, m
from ase.eos import EquationOfState
from ase.build import fcc111, fcc100, add_adsorbate
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
from ase.parallel import rank
from ase import Atoms
import sys
import os

homedir = os.path.expanduser('~')

experimental_lattice_parameter = 4.05

# # # Create Al bulk and initialize calculator parameters # # #
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
al_bulk = bulk('Al', 'fcc', a=experimental_lattice_parameter, cubic=False)

# # # Vary number of k-points to find a convergent number # # #
n_k_points = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

# # # Cutoff energy chosen somewhat arbitrarily, convergence checked in next program # # #
energy_cutoff = 200
energies = []
for n in n_k_points:
    if rank == 0:
        print 'Simulating ' + str(n) + 'kpoints...'
    mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
    calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                h=0.18,  # grid spacing, recommended value in this course
                xc='PBE',  # Exchange-correlation functional
                mixer=mixer,  # Mixer
                # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                kpts=(n, n, n),
                txt='simulate_bulk_Al_GPAW.txt')  # name of GPAW output text file

    al_bulk.set_calculator(calc)
    total_energy = al_bulk.get_potential_energy()
    energies.append(total_energy)

<<<<<<< HEAD
with open(homedir + '/TIF035/HA2/bulk/1_converge_kpoints_bulk.txt', 'w') as textfile:
=======
with open('/c3se/users/kasto/Hebbe/TIF035/HA2/bulk/1_converge_kpoints_bulk.txt', 'w') as textfile:
>>>>>>> 659d7c33b918c72c02078aa88bfcd2969a781cf7
    textfile.write('number of k points, bulk_energy\n')
    for i in range(len(n_k_points)):
        textfile.write(str(n_k_points[i]) + ',' + str(energies[i]) + '\n')
