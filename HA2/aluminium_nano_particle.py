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
import sys


class aluminium_nano_particle():
    def __init__(self):
        # # # Constants and parameters # # #
        self.a_al_experimental = 4.05  # Lattice constant for Al, experimentally determined
        self.E_ionization = 84.67567  # Ionization energy for hardest bound core electron
        self.surfaces = {}

    def save_nanoparticle():

    def load_nanoparticle():

    def create_bulk_aluminum(self, lattice_spacing):
        obj = bulk('Al', 'fcc', a=lattice_spacing, cubic=False)  # Create Al bulk
        self.bulk_Al = {'object': obj,
                        'a': lattice_spacing,
                        'volume': 0,
                        'energy': 0}

    def simulate_bulk_Al(self, E_cutoff, n_k_points):
        ''' Function that calculates volume and energy for bulk
            aluminium. The function uses GPAW with a plane wave basis set
            and the PBE exchange correlation. '''
        # # # Create Al bulk and initialize calculator parameters # # #

        mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
        calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                    h=0.18,  # grid spacing, recommended value in this course
                    xc='PBE',  # Exchange-correlation functional
                    mixer=mixer,
                    # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                    kpts=(n_k_points, n_k_points, n_k_points),
                    txt='simulate_bulk_Al_GPAW.txt')  # name of GPAW output text file

        self.bulk_Al['object'].set_calculator(calc)
        self.bulk_Al['volume'] = self.bulk_Al['object'].get_volume()
        self.bulk_Al['energy'] = self.bulk_Al['object'].get_potential_energy()

    def add_CO_adsorbate(self, surface, height_CO):
        ''' Add CO adsorbate on top of the Al bulk '''
        self.adsorbate = Atoms('CO')  # Create CO molecule object
        surface.add_adsorbate(slab=surf, adsorbate=self.adsorbate,
                              height=height_CO, position='ontop')
        # height above based on values for CO in ASE doc. Future: We could also
        # perform equilibrium scan by looping over various heights
        return 0  # Control value - 0 = code works, 1 = code does not

    def find_CO_opimal_distance(self, miller_index,):
        ''' Function that calculates the distance of the adsorbate which
        minimizes the energy '''
        surface = self.surfaces[miller_index]
        nbr_of_loops = 5
        theta = 1 / 4.0
        sigma_ads_vec = np.zeros(nbr_of_loops)
        h_vec = np.zeros(nbr_of_loops)
        for ind, h in enumerate(0.5, 5, nbr_of_loops):
            surface_tmp = self.surfaces[miller_index].copy()
            # sigma = self.simulate_surface_Al()
            self.add_CO_adsorbate(surface_tmp, h)
            self.adsorbate.set_cell([10, 10, 10])
            self.adsorbate.center()
            self.adsorbate.set_calculator(calc)
            energy_CO = self.adsorbate.get_potential_energy()  # Energy: CO
            surfEn_ads = surface_tmp.get_potential_energy()  # Energy: Surf with CO
            sigma = get_surface_sigma(miller_index)  # Sigma for bulk
            surfEn = get_surface_energy(miller_index)  # Energy: Surf w/o CO
            area = get_surface_area(miller_index)  # Surface area
            sigma_ads = sigma + theta * (surfEn_ads - surfEn - energy_CO) / area
            sigma_ads_vec[ind] = sigma_ads
            h_vec[ind] = h
        fig_1 = plt.figure()
        ax_potential = fig_1.add_subplot(111)
        ax_potential.plot(h_vec, sigma_ads_vec, label='Sigma with CO adsorbate')
        plt.savefig('eigAndEn.eps')
        plt.savefig('eigAndEn.png')
        return h_vec, sigma_ads_vec

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

    def simulate_surface_Al(self, miller_index):
        ''' Function that adds surfaces to the aluminum nano particle '''
        # Initialize new calculator that only considers k-space in xy-plane,
        # since we're only looking at the surface
        calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                    h=0.18,  # grid spacing
                    xc='PBE',  # XC-functional
                    mixer=mixer,
                    kpts=(k, k, 1),  # k-point grid
                    txt='surface_' + miller_index + '.txt')  # name of GPAW output text file

        self.get_surface_object(miller_index).set_calculator(calc)

        # Calc pot. energy of FCC
        surface_energy = self.surfaces[miller_index].get_potential_energy()

        area = get_surface_area(miller_index)
        # Calc. surf. energy per area (sigma) for FCC surface
        sigma = (1 / (2.0 * area)) * (surface_energy - N_x * N_y * E_bulk)

        self.surfaces[miller_index]['energy'] = surface_energy
        self.surfaces[miller_index]['sigma'] = sigma

    def calc_optimal_lattice_spacing(self, E_cutoff, n_k_points, n_points=7, q=0.02):
        ''' Function that calculates the lattice spacing that minimizes the
            potential energy for bulk aluminium. This is done by offsetting
            the experimentally derived lattice and then fitting the generated
            energies with GPAW's function EquationOfState.
            n_points: Number of lattice constants to loop over to find
                       equilibrium.
            q: Percentage difference from equalibrium for each point.'''
        energies = []
        volumes = []

        # # # Find lattice constant with lowest energy # # #
        cell_0 = al.cell  # Unit cell object of the Al bulk
        for eps in np.linspace(-q, q, n_points):
            self.bulk_Al['object'].cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell
            # Calculate the potential energy for the Al bulk
            simulate_bulk_Al(E_cutoff, n_k_points, lattice_spacing)
            energies.append(energy)
            volumes.append(volume)

            # Plot energies as a function of unit cell volume (directly related to latt. const.)
            eos = EquationOfState(volumes, energies)
            v0, E_bulk, B = eos.fit()
            eos.plot('Al_eos.png')
            # Latt. const. acc. to ASE doc., but why is this correct?
            a_calc = (4 * v0)**(1 / 3.0)

        if rank == 0:
            print 'Energies: ' + str(energies)
            print 'Volumes: ' + str(volumes)
            print 'a_calc: ' + str(a_calc)

        with open('calc_optimal_lattice_spacing.txt', 'w') as textfile:
            textfile.write('Optimal lattice spacing: ' + str(a_calc))
            textfile.write('Corresponding energy: ' + str(E_bulk))
            textfile.write('Corresponding bulk modulus: ' + str(B))
            textfile.write('Energies, Volumes')
            for i in range(len(volumes)):
                textfile.write(str(energies[i]) + ',' + str(volumes[i]))
        return a_calc


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
