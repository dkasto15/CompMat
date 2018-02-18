import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
from ase.build import fcc111, fcc100
import os

homedir = os.path.expanduser('~')

lattice_parameter = 4.04358928739
energy_bulk = -3.73707160907
n_k_points = 14
N_x = 1
N_y = 1
N_z = 7

slab111 = fcc111('Al', size=(N_x, N_y, N_z), a=lattice_parameter, vacuum=7.5)
cell111 = slab111.get_cell()  # Unit cell object of the Al FCC 111
area111 = np.linalg.norm(np.cross(cell111[0], cell111[1]))  # Calc. surface area

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

energy_cutoff = [50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600]
energies = []
sigmas = []

for energy in energy_cutoff:
    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n_k_points, n_k_points, 1),  # k-point grid
                txt='simulate_surface_Al_5_GPAW.txt')  # name of GPAW output text file

    slab111.set_calculator(calc)
    energy_slab = slab111.get_potential_energy()
    energies.append(energy_slab)
    sigmas.append((1 / (2.0 * area111)) * (energy_slab - N_z * energy_bulk))

with open(homedir + '/TIF035/HA2/surface/5_converge_cutoff_energy_surface.txt', 'w') as textfile:
    textfile.write('cutoff_energy, bulk_energy, surface_energy_density\n')
    for i in range(len(energy_cutoff)):
        textfile.write(str(energy_cutoff[i]) + ',' +
                       str(energies[i]) + ',' +
                       str(sigmas[i]) + '\n')
