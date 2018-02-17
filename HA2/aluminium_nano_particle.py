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
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import sys


class aluminium_nano_particle():
    def __init__(self):
        # # # Constants and parameters # # #
        self.experimental_lattice_parameter = 4.05  # Lattice constant for Al, experimentally determined
        self.E_ionization = 84.67567  # Ionization energy for hardest bound core electron
        self.surfaces = {}

    def create_bulk_Al(self, lattice_parameter):
        obj = bulk('Al', 'fcc', a=lattice_spacing, cubic=False)  # Create Al bulk
        self.bulk_Al = {'object': obj,
                        'lattice_parameter': lattice_spacing,
                        'volume': 0,
                        'energy': 0
                        'modulus': 0,
                        }

    def get_bulk_energy(self):
        return self.bulk_Al['energy']

    def get_bulk_lattice_parameter(self):
        return self.bulk_Al['lattice_parameter']

    def get_bulk_modulus(self):
        return self.bulk_Al['modulus']

    def get_bulk_cell(self):
        return self.bulk_Al['object'].cell

    def set_bulk_cell(self, cell):
        self.bulk_Al['object'].cell = cell

    def simulate_bulk_Al(self, energy_cutoff, n_k_points):
        ''' Function that calculates volume and energy for bulk
            aluminium. The function uses GPAW with a plane wave basis set
            and the PBE exchange correlation. '''
        # # # Create Al bulk and initialize calculator parameters # # #

        mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
        calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                    h=0.18,  # grid spacing, recommended value in this course
                    xc='PBE',  # Exchange-correlation functional
                    mixer=mixer,
                    # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                    kpts=(n_k_points, n_k_points, n_k_points),
                    txt='simulate_bulk_Al_GPAW.txt')  # name of GPAW output text file

        self.bulk_Al['object'].set_calculator(calc)
        self.bulk_Al['volume'] = self.bulk_Al['object'].get_volume()
        self.bulk_Al['energy'] = self.bulk_Al['object'].get_potential_energy()
        return self.bulk_Al['energy']

    def add_CO_adsorbate(self, miller_index, height_CO):
        ''' Add CO adsorbate on top of the Al bulk '''
        self.adsorbate = Atoms('CO')  # Create CO molecule object
        self.surfaces[miller_index]['object'].add_adsorbate(slab=surf, adsorbate=self.adsorbate,
                                                            height=height_CO, position='ontop')
        # height above based on values for CO in ASE doc. Future: We could also
        # perform equilibrium scan by looping over various heights
        return 0  # Control value - 0 = code works, 1 = code does not

    def simulate_CO(self, energy_cutoff, n_k_points):
        # # # Create Al bulk and initialize calculator parameters # # #

        mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
        calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                    h=0.18,  # grid spacing, recommended value in this course
                    xc='PBE',  # Exchange-correlation functional
                    mixer=mixer,
                    # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                    kpts=(n_k_points, n_k_points, n_k_points),
                    txt='simulate_CO_GPAW.txt')  # name of GPAW output text file

        self.adsorbate.set_calculator(calc)

        return self.adsorbate.get_potential_energy()

    def fixate_surface_atoms(self, miller_index):
        mask = [atom.symbol != 'Al' for atom in self.surfaces[miller_index]['object']]
        constraint = FixAtoms(mask=)
        self.surfaces[miller_index]['object'].set_constraint(constraint)

    def relax_adsorbate(self, miller_index):
        dyn = BFGS(atoms,
                   trajectory='relax_adsorbate_' + miller_index + '.traj',
                   logfile='relax_adsorbate_' + miller_index + '.traj')
        dyn.run(fmax=0.01)

    def create_surface_Al(self, miller_index, N_x, N_y, N_z, lattice_param):
        if miller_index == '111':
            obj = fcc111('Al', size=(N_x, N_y, N_z), a=lattice_param, vacuum=7.5)
            obj.set_calculator(calc)
            cell = obj.get_cell()  # Unit cell object of the Al FCC surface
            area = np.linalg.norm(np.cross(cell[0], cell[1]))  # Calc. surface area

            surface = {'object': obj,
                       'size': (N_x, N_y, N_z),
                       'a': lattice_param,
                       'area': area,
                       'energy': 0,
                       'sigma': 0}

        elif miller_index == '100':
            obj = fcc100('Al', size=(N_x, N_y, N_z), a=lattice_param, vacuum=7.5)
            obj.set_calculator(calc)
            cell = obj.get_cell()  # Unit cell object of the Al FCC surface
            area = np.linalg.norm(np.cross(cell[0], cell[1]))  # Calc. surface area

            surface = {'object': obj,
                       'size': (N_x, N_y, N_z),
                       'a': lattice_param,
                       'area': area,
                       'energy': 0,
                       'sigma': 0}
        else:
            return 1

        self.surface.center(axis=2)
        self.surfaces[miller_index].append(surface)
        return 0

    def get_surface_object(self, miller_index):
        return self.surfaces[miller_index]['object']

    def get_surface_area(self, miller_index):
        return self.surfaces[miller_index]['area']

    def get_surface_sigma(self, miller_index):
        return self.surfaces[miller_index]['sigma']

    def get_surface_energy(self, miller_index):
        return self.surfaces[miller_index]['energy']

    def simulate_surface_Al(self, miller_index, energy_cutoff, n_k_points):
        ''' Function that adds surfaces to the aluminum nano particle '''
        mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
        # Initialize new calculator that only considers k-space in xy-plane,
        # since we're only looking at the surface
        calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                    h=0.18,  # grid spacing
                    xc='PBE',  # XC-functional
                    mixer=mixer,
                    kpts=(n_k_points, n_k_points, 1),  # k-point grid
                    txt='simulate_surface_Al_' + miller_index + 'GPAW.txt')  # name of GPAW output text file

        self.get_surface_object(miller_index).set_calculator(calc)

        # Calc pot. energy of FCC
        surface_energy = self.surfaces[miller_index].get_potential_energy()
        self.surfaces[miller_index]['energy'] = surface_energy
        return surface_energy


with open('bulk_params.txt', 'w') as textfile:
    line = ''
    for el in E_cut:
        line = line + str(el) + ','
    textfile.write(line[:-1] + '\n')

    line = ''
    for el in k_points:
        line = line + str(el) + ','
    textfile.write(line[:-1] + '\n')

    for j in range(len(E_cut)):
        line = ''
        for i in range(len(k_points)):
            line = line + str(energies[i + j * len(k_points)]) + ','

        line = line[:-1]
        textfile.write(line + '\n')
