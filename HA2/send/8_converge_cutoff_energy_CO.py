#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW, Mixer, PW, PoissonSolver
# from gpaw.eigensolvers import Davidson
from ase.build import *
from ase import Atoms
from ase.parallel import rank
from ase.build import fcc111, fcc100
import os

homedir = os.path.expanduser('~')

d_CO = 1.128  # CO bondlength in [Å]
CO = Atoms('CO', positions=[(0., 0., 0.), (0., 0., d_CO)])
CO.set_cell([5, 5, 5])
CO.center()
#
# C = Atoms('C')
# C.set_cell([5, 5, 5])
# C.center()
#
# O = Atoms('O')
# O.set_cell([5, 5, 5])
# O.center()

mixer = Mixer(beta=0.25, nmaxold=3,	weight=1.0)  # Recommended values for small systems

energies = []
cutoff_energies = np.arange(300, 750, 50)

for cutoff_energy in cutoff_energies:
    if rank == 0:
        print 'Simulating ' + str(cutoff_energy) + ' cutoff...'

    calc = GPAW(mode=PW(cutoff_energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                txt='simulate_CO_8_GPAW.txt') # name of GPAW output text file

    # calc_C = GPAW(mode=PW(cutoff_energy),  # use the LCAO basis mode
    #             h=0.18,  # grid spacing
    #             xc='PBE',  # XC-functional
    #             mixer=mixer,
    #             eigensolver=Davidson(niter=5),
    #             txt='simulate_CO_8_GPAW_C.txt') # name of GPAW output text file

    CO.set_calculator(calc)
    energy_CO = CO.get_potential_energy()
    # C.set_calculator(calc_C)
    # energy_C = C.get_potential_energy()
    # O.set_calculator(calc)
    # energy_O = O.get_potential_energy()
    # energies.append(energy_CO - energy_C - energy_O)
    energies.append(energy_CO)

if rank == 0:
    plt.figure()
    plt.plot(cutoff_energies, energies, label='Energies')
    plt.xlabel('Cutoff energy [eV]')
    plt.ylabel('Energy CO [eV]')
    plt.savefig('cutoff_convergence_CO.png')
    plt.savefig('cutoff_convergence_CO.eps')

with open(homedir + '/TIF035/HA2/surface/8_converge_cutoff_CO.txt', 'w') as textfile:
    textfile.write('cutoff energy, CO energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' +
                       str(energies[i]) + '\n')
