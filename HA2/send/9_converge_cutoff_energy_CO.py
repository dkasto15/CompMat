import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
from ase.build import fcc111, fcc100
import os

homedir = os.path.expanduser('~')

CO = Atoms('CO')

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

n_k_points = 16
energy_cutoff = 300
energies = []

n_k_points = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
for n in n_k_points:
    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n, n, 1),  # k-point grid
                txt='simulate_CO_8_GPAW.txt')  # name of GPAW output text file

    CO.set_calculator(calc)
    energy_pot = CO.get_potential_energy()
    energies.append(energy_pot)

with open(homedir + '/TIF035/HA2/surface/8_converge_kpoints_CO.txt', 'w') as textfile:
    textfile.write('atom depth, bulk_energy\n')
    for i in range(len(N_z_vec)):
        textfile.write(str(N_z_vec[i]) + ',' +
                       str(energies[i]) + '\n')