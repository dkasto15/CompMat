import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank
from ase.eos import EquationOfState
import os

homedir = os.path.expanduser('~')

n_k_points = 16
energy_cutoff = 400

q = 0.02 # Percental difference from experimental lattice parameter for each point
n_points = 7 # n_points: Number of lattice constants to loop over to find equilibrium.

experimental_lattice_parameter = 4.05

# # # Create Al bulk and initialize calculator parameters # # #
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
al_bulk = bulk('Al', 'fcc', a=experimental_lattice_parameter, cubic=False)

calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
            h=0.18,  # grid spacing, recommended value in this course
            xc='PBE',  # Exchange-correlation functional
            mixer=mixer,  # Mixer
            # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
            kpts=(n_k_points, n_k_points, n_k_points),
            txt='simulate_bulk_Al_GPAW_2.txt')  # name of GPAW output text file
al_bulk.set_calculator(calc)

energies = []
volumes = []

# # # Find lattice constant with lowest energy # # #
cell_0 = al_bulk.cell  # Unit cell object of the Al bulk

for eps in np.linspace(-q, q, n_points):

    al_bulk.cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell

    # Calculate the potential energy for the Al bulk
    total_energy = al_bulk.get_potential_energy()
    energies.append(total_energy)

    volume = al_bulk.get_volume()
    volumes.append(volume)

# Plot energies as a function of unit cell volume (directly related to latt. const.)
eos = EquationOfState(volumes, energies)
v0, E_bulk, B = eos.fit()
eos.plot('Al_eos.png')

# Latt. const. acc. to ASE doc., but why is this correct?
a_calc = (4 * v0)**(1 / 3.0)


with open(homedir + '/TIF035/HA3/bulk/3_calc_optimal_lattice_spacing.txt', 'w') as textfile:
    textfile.write('Optimal lattice spacing: ' + str(a_calc))
    textfile.write('Corresponding energy: ' + str(E_bulk))
    textfile.write('Corresponding bulk modulus: ' + str(B))
    textfile.write('Energies, Volumes')
    for i in range(len(volumes)):
        textfile.write(str(energies[i]) + ',' + str(volumes[i]))
