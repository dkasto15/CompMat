
#from ase.io import read
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np

# input files
filename1 = '7_calc_adsorbtion_energy.txt'
filename2 = '7_surface_sigma.txt'

E_CO = -13.9
theta = 1/4

# import and manage data
# sigma_ads = []

with open('HA2/' + filename1, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        if line[0] == '111':
            sigma_ads_111 = float(line[3])
            area_111 = float(line[1])
        if line[0] == '100':
            sigma_ads_100 = float(line[3])
            area_100 = float(line[1])

with open('HA2/' + filename2, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        if line[0] == '111':
            sigma_111 = float(line[4])
        if line[0] == '100':
            sigma_100 = float(line[4])

sigma_int_111 = sigma_111 + theta * (sigma_ads_111 - sigma_111 - E_CO / area_111)
sigma_int_100 = sigma_100 + theta * (sigma_ads_100 - sigma_100 - E_CO / area_100)

print(sigma_int_111)
print(sigma_int_100)


''' Case: Just Al '''
Al = wulff_construction('Al',
                        surfaces=[(1, 0, 0),
                                  (1, 1, 1)],
                        energies=[sigma_100, sigma_111],
                        size=10000,
                        structure='fcc',
                        rounding='below')  # Round to closest structure below size
Al.center(vacuum=10)
view(Al)

''' Case: Al + adsorbed CO '''
Al_ads = wulff_construction('Al',
                            surfaces=[(1, 0, 0),
                                      (1, 1, 1)],
                            energies=[sigma_int_100, sigma_int_111],
                            size=10000,
                            structure='fcc',
                            rounding='below')  # Round to closest structure below size
Al_ads.center(vacuum=10)
view(Al_ads)
