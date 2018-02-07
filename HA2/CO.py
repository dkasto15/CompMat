
#from ase.io import read
import numpy as np
#from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.io import read
from ase.units import kJ, J, m
from ase.eos import EquationOfState
from ase.build import fcc111, fcc100, add_adsorbate
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
from ase.parallel import rank
from ase import Atoms

N_x = 2
N_y = 2
N_z = 20

a_al = 4.05

surface111 = fcc111('Al', size=(N_x, N_y, N_z), a=a_al, vacuum=7.5)
surface100 = fcc100('Al', size=(N_x, N_y, N_z), a=a_al, vacuum=7.5)
surface111.center(axis=2)
surface100.center(axis=2)
CO = Atoms('CO')
add_adsorbate(slab=surface111, adsorbate=CO, height=4.5, position='ontop')


# view(surface100)
view(surface111)

# Al = wulff_construction('Al',
#                            surfaces=[(1, 0, 0),
#                                      (1, 1, 1)],
#                            energies=[sigma100, sigma111],
#                            size=10000,
#                            structure='fcc',
#                            rounding='below')  # What does this one do?
# Al.center(vacuum=10)
# view(Al)

# Al_ads = wulff_construction('Al',
#                            surfaces=[(1, 0, 0),
#                                      (1, 1, 1)],
#                            energies=[sigma100_ads, sigma111_ads],
#                            size=10000,
#                            structure='fcc',
#                            rounding='below')  # What does this one do?
# Al_ads.center(vacuum=10)
# view(Al_ads)
