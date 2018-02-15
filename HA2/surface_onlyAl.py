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

if('layers' not in sys.argv):
    sys.exit('\"layers\" not in argument variable')
n_layers = sys.argv.index('layers')

if('cutoff' not in sys.argv):
    sys.exit('\"cutoff\" not in argument variable')
n_cut = sys.argv.index('cutoff')

if('kspace' not in sys.argv):
    sys.exit('\"kspace\" not in argument variable')
n_k = sys.argv.index('kspace')


try:
    nbrOfLayers = float(sys.argv[n_layers + 1])
    E_cut = [float(el) for el in sys.argv[n_cut + 1:n_k]]
    k_points = [float(el) for el in sys.argv[n_k + 1:]]
except Exception as e:
    print e
    sys.exit('Format: python ' + sys.argv[0] + 'layers + n'
             'cutoff E1 E2 ... ' + 'kspace k1 k2 ...')

# # # Create surface slab with one atom and nbrOfLayers layers # # #
N_x = 1
N_y = 1
N_z = nbrOfLayers

surface111 = fcc111('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
surface100 = fcc100('Al', size=(N_x, N_y, N_z), a=a_calc, vacuum=7.5)
surface111.center(axis=2)
surface100.center(axis=2)

for energy in E_cut:
    for k in k_points:
        # Initialize new calculator that only considers k-space in xy-plane,
        # since we're only looking at the surface
        calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                    h=0.18,  # grid spacing
                    xc='PBE',  # XC-functional
                    mixer=mixer,
                    kpts=(k, k, 1),  # k-point grid
                    txt='surface_onlyAl_GPAW.txt')  # name of GPAW output text file

        surface111.set_calculator(calc)
        surface100.set_calculator(calc)

        cell111 = surface111.get_cell()  # Unit cell object of the Al FCC 111
        area111 = np.linalg.norm(np.cross(cell111[0], cell111[1]))  # Calc. surface area
        surfEn111 = surface111.get_potential_energy()  # Calc pot. energy of FCC 111
        cell100 = surface100.get_cell()  # Unit cell object of the Al FCC 100
        area100 = np.linalg.norm(np.cross(cell100[0], cell100[1]))  # Calc. surface area
        surfEn100 = surface100.get_potential_energy()  # Calc pot. energy of FCC 100

        # Calc. surf. energy per area (sigma) for FCC 111 and 100
        sigma111 = (1 / (2.0 * area111)) * (surfEn111 - N_x * N_y * E_bulk)
        sigma100 = (1 / (2.0 * area100)) * (surfEn100 - N_x * N_y * E_bulk)


with open('surface_onlyAl_params.txt', 'w') as textfile:
    textfile.write(str(k_points) + '\n')
    textfile.write(str(E_cut) + '\n')

    for j in range(len(E_cut)):
        line = ''
        for i in range(len(k_points)):
            line = line + str(energies[i + j * len(k_points)]) + ','

        line = line[:-1]
        textfile.write(line + '\n')
