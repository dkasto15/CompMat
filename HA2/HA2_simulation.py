#!/usr/bin/env python

# # imports # #
import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import bulk
from ase.io import read
from ase.units import kJ
from ase.utils.eos import EquationOfState
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view

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
for eps in np.linspace(-0.02,	0.02,	N_lattice_spacings):
    al.cell = (1 + eps) * cell_0
    al.get_potential_energy()
    print al.get_potential_energy()
