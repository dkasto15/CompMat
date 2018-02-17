from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
import os

homedir = os.path.expanduser('~')

n_k_points = 16
experimental_lattice_parameter = 4.05

# # # Create Al bulk and initialize calculator parameters # # #
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
al_bulk = bulk('Al', 'fcc', a=experimental_lattice_parameter, cubic=False)

energies = []
cutoff_energies = [20, 50, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600]
for cutoff_energy in cutoff_energies:
    mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
    calc = GPAW(mode=PW(cutoff_energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing, recommended value in this course
                xc='PBE',  # Exchange-correlation functional
                mixer=mixer,  # Mixer
                # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                kpts=(n_k_points, n_k_points, n_k_points),
                txt='simulate_bulk_Al_GPAW_2.txt')  # name of GPAW output text file

    al_bulk.set_calculator(calc)
    total_energy = al_bulk.get_potential_energy()
    energies.append(total_energy)

with open(homedir + '/TIF035/HA2/bulk/2_converge_cutoff_energy_bulk.txt', 'w') as textfile:
    textfile.write('cutoff_energy, bulk_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' + str(energies[i]) + '\n')
