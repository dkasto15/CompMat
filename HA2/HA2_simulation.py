#!/usr/bin/env python
# coding=utf-8

# # imports # #
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

#from ase.lattice.surface import *

a_al = 4.05
al = bulk('Al', 'fcc', a=a_al, cubic=False)

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)
N_lattice_spacings = 7
E_al = 84.67567  # ionization energy for hardest bound core electron
E_cut = [50, E_al, 100, 200, 300, 400, 500, 600, 700, 800]  # cut-off energy

for energy in [500]:
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(12, 12, 12),  # k-point grid
                txt='out.txt')  # name of GPAW output text file

    al.set_calculator(calc)
    cell_0 = al.cell
    for eps in np.linspace(-0.02, 0.02,	N_lattice_spacings):
        al.cell = (1 + eps) * cell_0
        al.get_potential_energy()  # ghj

    confs = read('out.txt@0:7')  # read 7 conﬁgurations

    # Extract volumes and energies:
    volumes = [atoms.get_volume() for atoms in confs]
    energies = [atoms.get_potential_energy() for atoms in confs]
    # if rank == 0:
    #     print energies, shape(energies)
    eos = EquationOfState(volumes, energies)
    v0, E_bulk, B = eos.fit()
    eos.plot('Al_eos.png')
    a_calc = (4 * v0)**(1 / 3.0)  # Is this correct?

    N_x = 1
    N_y = 1
    N_z = 10

    surface111 = fcc111('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
    surface100 = fcc100('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
    surface111.center(axis=2)
    surface100.center(axis=2)

    calc2 = GPAW(mode=PW(energy),  # use the LCAO basis mode
                 h=0.18,  # grid spacing
                 xc='PBE',  # XC-functional
                 mixer=mixer,
                 kpts=(12, 12, 1),  # k-point grid
                 txt='out2.txt')  # name of GPAW output text file

    surface111.set_calculator(calc2)
    surface100.set_calculator(calc2)

    cell111 = surface111.get_cell()
    area111 = np.linalg.norm(np.cross(cell111[0], cell111[1]))
    surfEn111 = surface111.get_potential_energy()

    cell100 = surface100.get_cell()
    area100 = np.linalg.norm(np.cross(cell100[0], cell100[1]))
    surfEn100 = surface100.get_potential_energy()

    sigma111 = (1 / (2.0 * area111)) * (surfEn111 - N_x * N_y * E_bulk)
    sigma100 = (1 / (2.0 * area100)) * (surfEn100 - N_x * N_y * E_bulk)

    file = open('sigmas.txt', 'w')
    file.write(str(sigma111) + '\t' + str(sigma100))
    file.close()

    # # # Add adsorbate # # #
    d_CO = 1.128  # CO bondlength in [Å]
    CO = Atoms('CO', [(0, 0, 0), (0, 0, d_CO)])
    CO
    add_adsorbate(slab=surface111, adsorbate=CO, height=4.5, position='ontop')
    add_adsorbate(slab=surface100, adsorbate=CO, height=4.5, position='ontop')

    cell111 = surface111.get_cell()
    area111 = np.linalg.norm(np.cross(cell111[0], cell111[1]))
    surfEn111 = surface111.get_potential_energy()

    cell100 = surface100.get_cell()
    area100 = np.linalg.norm(np.cross(cell100[0], cell100[1]))
    surfEn100 = surface100.get_potential_energy()

    sigma111_ads = (1 / (2.0 * area111)) * (surfEn111 - N_x * N_y * E_bulk)
    sigma100_ads = (1 / (2.0 * area100)) * (surfEn100 - N_x * N_y * E_bulk)

    file = open('sigmas_ads.txt', 'w')
    file.write(str(sigma111_ads) + '\t' + str(sigma100_ads))
    file.close()

    if rank == 0:
        a = 0
        print area111, surfEn111, E_bulk, sigma111
