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

    def simulate_bulk_Al(self, E_cutoff, n_k_points, lattice_spacing):
        ''' Function that calculates volume and energy for bulk
            aluminium. The function uses GPAW with a plane wave basis set
            and the PBE exchange correlation. '''

        # # # Create Al bulk and initialize calculator parameters # # #
        al = bulk('Al', 'fcc', a=lattice_spacing, cubic=False)  # Create Al bulk
        mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems

        calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                    h=0.18,  # grid spacing, recommended value in this course
                    xc='PBE',  # Exchange-correlation functional
                    mixer=mixer,
                    # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                    kpts=(n_k_points, n_k_points, n_k_points),
                    txt='simulate_bulk_Al_GPAW.txt')  # name of GPAW output text file
        al.set_calculator(calc)
        return (al.get_volume(), al.get_potential_energy())

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
                al.cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell
                # Calculate the potential energy for the Al bulk
                (energy, volume) = simulate_bulk_Al(E_cutoff, n_k_points, lattice_spacing):
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
                textfile.write('Correspongng bulk modulus: ' + str(B))
                textfile.write('Energies, Volumes')
                for i in range(len(volumes)):
                    textfile.write(str(energies[i]) + ',' + str(volumes[i]))
            return a_calc

        def simulate_surface_Al(surface, k):


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
