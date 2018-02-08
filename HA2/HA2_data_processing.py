
#from ase.io import read
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np

# input files
filename1 = 'sigma_Al.txt'
filename2 = 'sigma_ads.txt'

# import and manage data
sigma_Al = np.array(np.loadtxt(filename1))
sigma111_Al = sigma_Al[0]
sigma100_Al = sigma_Al[1]
sigmas_ads = np.array(np.loadtxt(filename2))
sigma111_ads = sigmas_ads[0]
sigma100_ads = sigmas_ads[1]

Al = wulff_construction('Al',
                        surfaces=[(1, 0, 0),
                                  (1, 1, 1)],
                        energies=[sigma100_Al, sigma111_Al],
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
