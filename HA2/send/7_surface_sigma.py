import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
from ase.build import fcc111, fcc100
import os

homedir = os.path.expanduser('~')

lattice_parameter = 4.04358928739
energy_bulk = -3.73707160907

N_x = 1
N_y = 1
N_z = 14

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for clusters

n_k_points = 16
energy_cutoff = 350

energies = []
sigmas = []
areas = []

surfaces = []
surfaces.append(fcc111('Al', size=(N_x, N_y, N_z), a=lattice_parameter, vacuum=7.5))
surfaces.append(fcc100('Al', size=(N_x, N_y, N_z), a=lattice_parameter, vacuum=7.5))
miller_indices = ['111', '100']

for i, slab in enumerate(surfaces):
    if rank == 0:
        print 'Simulating ' + miller_indices[i] + '...'
    slab.center(axis=2)
    cell = slab.get_cell()  # Unit cell object
    area = np.linalg.norm(np.cross(cell[0], cell[1]))  # Calc. surface area
    areas.append(area)
    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n_k_points, n_k_points, 1),  # k-point grid
                txt='surface_sigma_7_GPAW_' + str(miller_indices[i]) + '.txt')  # name of GPAW output text file
    slab.set_calculator(calc)
    energy_slab = slab.get_potential_energy()
    energies.append(energy_slab)
    sigmas.append((1 / (2.0 * area)) * (energy_slab - N_z * energy_bulk))



with open(homedir + '/TIF035/HA2/surface/7_surface_sigma.txt', 'w') as textfile:
    textfile.write('miller index, atom depth, area, bulk energy, surface energy density\n')
    for i in range(len(surfaces)):
        textfile.write(str(miller_indices[i]) + ',' +
                       str(N_z) + ',' +
                       str(areas[i]) + ',' +
                       str(energies[i]) + ',' +
                       str(sigmas[i]) + '\n')
