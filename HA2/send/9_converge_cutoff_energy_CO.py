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

CO = Atoms('CO')

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

n_k_points = 16

energies = []
cutoff_energies = [50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600]

for cutoff_energy in cutoff_energies:
    if rank == 0:
        print 'Simulating ' + str(cutoff_energy) + ' cutoff...'

    calc = GPAW(mode=PW(cutoff_energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n_k_points, n_k_points, 1),  # k-point grid
                txt='simulate_CO_8_GPAW.txt')  # name of GPAW output text file

    CO.set_calculator(calc)
    energy_pot = CO.get_potential_energy()
    energies.append(energy_pot)

    fig_2 = plt.figure()
    fig_2.plot(cutoff_energies, energies, label='Energies')
    fig_2.set_xlabel('Cutoff energy [a.u.]')
    fig_2.set_ylabel('Energy [atomic units]')
    fig_2.savefig('cutoff_convergence_CO.png')
    fig_2.savefig('cutoff_convergence_CO.eps')

with open(homedir + '/TIF035/HA2/surface/9_converge_cutoff_CO.txt', 'w') as textfile:
    textfile.write('cutoff energy, CO_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' +
                       str(energies[i]) + '\n')
