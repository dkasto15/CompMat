
#from ase.io import read
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np

# input files
filename1 = 'sigmas.txt'
filename2 = 'sigmas_ads.txt'

# import and manage data
sigmas = np.array(np.loadtxt(filename1))
sigma111 = sigmas[0]
sigma100 = sigmas[1]
sigmas_ads = np.array(np.loadtxt(filename2))
sigma111_ads = sigmas_ads[0]
sigma100_ads = sigmas_ads[1]

print(sigma111, sigma100)
print(sigma111_ads, sigma100_ads)

Al = wulff_construction('Al',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1)],
                           energies=[sigma100, sigma111],
                           size=10000,
                           structure='fcc',
                           rounding='below')  # What does this one do?
Al.center(vacuum=10)
view(Al)

Al_ads = wulff_construction('Al',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1)],
                           energies=[sigma100_ads, sigma111_ads],
                           size=10000,
                           structure='fcc',
                           rounding='below')  # What does this one do?
Al_ads.center(vacuum=10)
view(Al_ads)
