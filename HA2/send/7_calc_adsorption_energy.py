#!/usr/bin/env python
# coding=utf-8
import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
from ase.build import fcc111, fcc100, add_adsorbate
from ase import Atoms
from ase.optimize import BFGS, QuasiNewton
from ase.constraints import FixAtoms
from ase.io import write
from gpaw.poisson import PoissonSolver
import os

homedir = os.path.expanduser('~')
file_simulation_parameters = open('7_simulation_params.txt', 'w')
lattice_parameter = 4.04358928739
energy_bulk = -3.73707160907

N_x = 1
N_y = 1
N_z = 14
n_k_points = 16
energy_cutoff = 350

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

surfaces = []
surfaces.append(fcc111('Al', size=(N_x, N_y, N_z), a=lattice_parameter, vacuum=15))
surfaces.append(fcc100('Al', size=(N_x, N_y, N_z), a=lattice_parameter, vacuum=15))
miller_indices = ['111', '100']

d_CO = 1.128  # CO bondlength in [Å]
CO_adsorbate = Atoms('CO', positions=[(0., 0., 0.), (0., 0., d_CO)])

energies = []
sigmas = []
areas = []

for i, slab in enumerate(surfaces):
    cell = slab.get_cell()  # Unit cell object of the Al FCC 111
    area = np.linalg.norm(np.cross(cell[0], cell[1]))  # Calc. surface area
    areas.append(area)

    slab.center(axis=2)

    add_adsorbate(slab=slab, adsorbate=CO_adsorbate,
                  height=2, position='ontop')

    write('slab' + miller_indices[i] + '.png', slab, rotation='10z,-80x')
    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(n_k_points, n_k_points, 1),  # k-point grid
                txt='simulate_surface_Al_7_GPAW.txt')  # name of GPAW output text file

    slab.set_calculator(calc)

    mask = [atom.symbol == 'Al' for atom in slab]

    file_simulation_parameters.write('mask' + '\n')
    for el in mask:
        file_simulation_parameters.write(str(el) + '\n')
    constraint = FixAtoms(mask=mask)
    slab.set_constraint(constraint)

    dyn = BFGS(slab,
               trajectory='relax_adsorbate_' + miller_indices[i] + '.traj',
               logfile='relax_adsorbate_' + miller_indices[i] + '.qn')
    dyn.run(fmax=0.01)

    energy_slab = slab.get_potential_energy()
    energies.append(energy_slab)
    sigmas.append((1 / (2.0 * area)) * (energy_slab - N_z * energy_bulk))

file_simulation_parameters.close()
with open(homedir + '/TIF035/HA2/surface/7_calc_adsorbtion_energy.txt', 'w') as textfile:
    textfile.write('miller_index, area, bulk_energy, surface_energy_density\n')
    for i in range(len(surfaces)):
        textfile.write(str(miller_indices[i]) + ',' +
                       str(areas[i]) + ',' +
                       str(energies[i]) + ',' +
                       str(sigmas[i]) + '\n')
