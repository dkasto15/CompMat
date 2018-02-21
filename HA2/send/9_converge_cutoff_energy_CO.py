#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW, Mixer, PW, PoissonSolver
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

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

energies = []
cutoff_energies = [100, 150, 200, 250, 300, 350, 400, 450, 500, 600]

for cutoff_energy in cutoff_energies:
    if rank == 0:
        print 'Simulating ' + str(cutoff_energy) + ' cutoff...'

    calc = GPAW(mode=PW(cutoff_energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                txt='simulate_CO_9_GPAW.txt') # name of GPAW output text file

    CO.set_calculator(calc)
    energy_pot = CO.get_potential_energy()
    energies.append(energy_pot)

if rank == 0:
    plt.figure()
    plt.plot(cutoff_energies, energies, label='Energies')
    plt.xlabel('Cutoff energy [eV]')
    plt.ylabel('Energy [atomic units]')
    plt.savefig('cutoff_convergence_CO.png')
    plt.savefig('cutoff_convergence_CO.eps')

with open(homedir + '/TIF035/HA2/surface/9_converge_cutoff_CO.txt', 'w') as textfile:
    textfile.write('cutoff energy, CO_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' +
                       str(energies[i]) + '\n')
