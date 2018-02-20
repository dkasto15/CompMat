#!/usr/bin/env python
# coding=utf-8

import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase import Atoms
from ase.parallel import rank
from ase.build import fcc111, fcc100
import os

homedir = os.path.expanduser('~')

d_CO = 1.128  # CO bondlength in [Å]
CO = Atoms('CO', positions=[(0., 0., 0.), (0., 0., d_CO)])
CO.set_cell([10, 10, 10])
CO.center()

mixer = Mixer(beta=0.25, nmaxold=3, weight=1.0)  # Recommended values for small systems

energy_cutoff = 350
energies = []

n_k_points = range(4, 21)
for n in n_k_points:
    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    if rank == 0:
        print 'Simulating ' + str(n) + ' k_point...'
    calc = GPAW(mode='lcao',  # use the LCAO basis mode
                h=0.18,  # grid spacing
                basis='dzp',
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n, n, n),  # k-point grid
                txt='simulate_CO_8_GPAW.txt')  # name of GPAW output text file

    CO.set_calculator(calc)
    energy_pot = CO.get_potential_energy()
    energies.append(energy_pot)

    fig_2 = plt.figure()
    fig_2.plot(n_k_points, energies, label='Energies')
    fig_2.set_xlabel('Number of k points [a.u.]')
    fig_2.set_ylabel('Energy [atomic units]')
    fig_2.savefig('k_point_convergence.png')
    fig_2.savefig('k_point_convergence.eps')

with open(homedir + '/TIF035/HA2/surface/8_converge_kpoints_CO.txt', 'w') as textfile:
    textfile.write('Number of k-points, CO energy\n')
    for i in range(len(n_k_points)):
        textfile.write(str(n_k_points[i]) + ',' +
                       str(energies[i]) + '\n')
