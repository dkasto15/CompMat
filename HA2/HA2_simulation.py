#!/usr/bin/env python
#coding=utf-8

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

# # # Constants and parameters # # #
a_al = 4.05 # Lattice constant for Al, experimentally determined
N_lattice_spacings = 7 # Number of lattice constants to loop over to find equilibrium
E_al = 84.67567  # Ionization energy for hardest bound core electron

# # # Create Al bulk and initialize calculator parameters # # #
al = bulk('Al', 'fcc', a=a_al, cubic=False) # Create Al bulk
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0) # Recommended values for small systems
E_cut = [50, E_al, 100, 200, 300, 400, 500]  # cut-off energy

for energy in [500]: # Change to E_cut to loop and check convergence
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing, recommended value in this course
                xc='PBE',  # Exchange-correlation functional
                mixer=mixer,
                kpts=(12, 12, 12),  # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                txt='out.txt')  # name of GPAW output text file
    al.set_calculator(calc)

    # # # Find lattice constant with lowest energy # # #
    cell_0 = al.cell # Unit cell object of the Al bulk
    for eps in np.linspace(-0.02, 0.02,	N_lattice_spacings):
        al.cell = (1 + eps) * cell_0 # Adjust lattice constant of unit cell
        al.get_potential_energy() # Calculate the potential energy for the Al bulk

    confs = read('out.txt@0:'+str(N_lattice_spacings))  # Read the conﬁgurations

    # Extract volumes and energies:
    volumes = [atoms.get_volume() for atoms in confs]
    energies = [atoms.get_potential_energy() for atoms in confs]
    # if rank == 0:
    #     print energies, shape(energies)

    # Plot energies as a function of unit cell volume (directly related to latt. const.)
    eos = EquationOfState(volumes, energies)
    v0, E_bulk, B = eos.fit()
    eos.plot('Al_eos.png')
    a_calc = (4 * v0)**(1 / 3.0) # Latt. const. acc. to ASE doc., but why is this correct?

    # # # Create Al FCC structures (111 and 100) with 2x2 atoms on the
    #     surface and 10 layers # # #
    N_x = 2
    N_y = 2
    N_z = 10

    surface111 = fcc111('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
    surface100 = fcc100('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
    surface111.center(axis=2)
    surface100.center(axis=2)

    # Initialize new calculator that only considers k-space in xy-plane,
    # since we're only looking at the surface
    calc2 = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(12, 12, 1),  # k-point grid
                txt='out2.txt')  # name of GPAW output text file

    surface111.set_calculator(calc2)
    surface100.set_calculator(calc2)

    cell111 = surface111.get_cell() # Unit cell object of the Al FCC 111
    area111 = np.linalg.norm(np.cross(cell111[0], cell111[1])) # Calc. surface area
    surfEn111 = surface111.get_potential_energy() # Calc pot. energy of FCC 111
    cell100 = surface100.get_cell() # Unit cell object of the Al FCC 100
    area100 = np.linalg.norm(np.cross(cell100[0], cell100[1])) # Calc. surface area
    surfEn100 = surface100.get_potential_energy() # Calc pot. energy of FCC 100

    # Calc. surf. energy per area (sigma) for FCC 111 and 100
    sigma111 = (1 / (2.0 * area111)) * (surfEn111 - N_x * N_y * E_bulk)
    sigma100 = (1 / (2.0 * area100)) * (surfEn100 - N_x * N_y * E_bulk)

    # Save sigmas for 111 and 100 to file
    file = open('sigmas.txt','w')
    file.write(str(sigma111)+'\t'+str(sigma100))
    file.close()

    # # # Add CO adsorbate to Al surface # # #
    d_CO = 1.128  # CO bondlength in [Å]
    CO = Atoms('CO') # Create CO molecule object
    add_adsorbate(slab=surface111, adsorbate=CO, height=4.5, position='ontop')
    add_adsorbate(slab=surface100, adsorbate=CO, height=4.5, position='ontop')
    # height above based on values in ASE doc. Future: We could also perform equilibrium
    # scan by looping over various heights


    cell111 = surface111.get_cell()
    area111 = np.linalg.norm(np.cross(cell111[0], cell111[1]))
    surfEn111 = surface111.get_potential_energy()
    cell100 = surface100.get_cell()
    area100 = np.linalg.norm(np.cross(cell100[0], cell100[1]))
    surfEn100 = surface100.get_potential_energy()

    sigma111_ads = (1 / (2.0 * area111)) * (surfEn111 - N_x * N_y * E_bulk)
    sigma100_ads = (1 / (2.0 * area100)) * (surfEn100 - N_x * N_y * E_bulk)

    file = open('sigmas_ads.txt','w')
    file.write(str(sigma111_ads)+'\t'+str(sigma100_ads))
    file.close()

    if rank == 0:
        a = 0
        print area111, surfEn111, E_bulk, sigma111
