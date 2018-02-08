#!/usr/bin/env python
# coding=utf-8

# # # imports # # #
import numpy as np
#from gpaw import GPAW, Mixer, PW
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

# # # Constants and parameters # # #
a_al = 4.05  # Lattice constant for Al, experimentally determined
N_lattice_spacings = 7  # Number of lattice constants to loop over to find equilibrium
E_al = 84.67567  # Ionization energy for hardest bound core electron

# # # Create Al bulk and initialize calculator parameters # # #
al = bulk('Al', 'fcc', a=a_al, cubic=False)  # Create Al bulk
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
volumes = []
energies = []
if('loop_lattice_param' in sys.argv):
    loop_lattice_param = True
else:
    loop_lattice_param = False

if('cutoff' not in sys.argv):
    sys.exit('\"cutoff\" not in argument variable')
n_cut = sys.argv.index('cutoff')

if('kspace' not in sys.argv):
    sys.exit('\"kspace\" not in argument variable')
n_k = sys.argv.index('kspace')


try:
    E_cut = [float(el) for el in sys.argv[n_cut + 1:n_k]]
    k_points = [float(el) for el in sys.argv[n_k + 1:]]
except Exception as e:
    print e
    sys.exit('Format: python loop_lattice_param' + sys.argv[0] +
             'cutoff E1 E2 ... ' + 'kspace k1 k2 ...')


for energy in E_cut:
    for k in k_points:
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing, recommended value in this course
                xc='PBE',  # Exchange-correlation functional
                mixer=mixer,
                kpts=(k, k, k),  # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                txt='out.txt')  # name of GPAW output text file
    if loop_lattice_param:
        # # # Find lattice constant with lowest energy # # #
        cell_0 = al.cell  # Unit cell object of the Al bulk
        for eps in np.linspace(-0.02, 0.02,	N_lattice_spacings):
            al.cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell
            # Calculate the potential energy for the Al bulk
            energies.append(al.get_potential_energy())
            volumes.append(al.get_volume())

        #
        # confs = read('out.txt@0:' + str(N_lattice_spacings))  # Read the conﬁgurations
        # # Extract volumes and energies:
        # volumes = [atoms.get_volume() for atoms in confs]
        # energies = [atoms.get_potential_energy() for atoms in confs]
        # # if rank == 0:
        # #     print energies, shape(energies)
        #
        if rank == 0:
            print 'Energies: ' + str(energies)
            print 'Volumes: ' + str(volumes)
            # Plot energies as a function of unit cell volume (directly related to latt. const.)
            eos = EquationOfState(volumes, energies)
            v0, E_bulk, B = eos.fit()
            eos.plot('Al_eos.png')
            a_calc = (4 * v0)**(1 / 3.0)  # Latt. const. acc. to ASE doc., but why is this correct?

            print 'a_calc: ' + str(a_calc)

    else:
        al.get_potential_energy()
