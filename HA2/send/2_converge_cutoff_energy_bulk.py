from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank

try:
    with open('~/TIF035/HA2/bulk/converged_kpoint.txt', 'r') as textfile:
        n_k_points = int(next(textfile))
except Exception as e:
    print('converged_kpoint.txt does not exist')
    exit(1)

experimental_lattice_parameter = 4.05

# # # Create Al bulk and initialize calculator parameters # # #
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
al_bulk = bulk('Al', 'fcc', a=experimental_lattice_parameter, cubic=False)

energies = []
cutoff_energies = [50, 100, 200, 300, 400, 500]
for cutoff_energy in cutoff_energies:
    particle.simulate_bulk_Al(cutoff_energy, n_k_points)
    energies.append(particle.get_bulk_energy())

with open('~/TIF035/HA2/bulk/2_converge_cutoff_energy_bulk.txt', 'w') as textfile:
    texfile.write('cutoff_energy, bulk_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' + str(energies[i]))
