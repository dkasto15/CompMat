
#from ase.io import read
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np

# input file
filename1 = 'sigmas.txt'

# import and manage data
sigmas = np.array(np.loadtxt(filename1))
sigma111 = sigmas[0]
sigma100 = sigmas[1]

print(sigma111, sigma100)

atoms = wulff_construction('Al',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1)],
                           energies=[sigma100, sigma111],
                           size=10000,
                           structure='fcc',
                           rounding='below')  # What does this one do?
atoms.center(vacuum=10)
view(atoms)
