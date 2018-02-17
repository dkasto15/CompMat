from gpaw import GPAW, Mixer, PW
from ase.build import *
from ase.parallel import rank

experimental_lattice_parameter = 4.05

# # # Create Al bulk and initialize calculator parameters # # #
mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
al_bulk = bulk('Al', 'fcc', a=experimental_lattice_parameter, cubic=False)

# # # Vary number of k-points to find a convergent number # # #
n_k_points = [4, 5, 6, 7, 8, 9, 10]

# # # Cutoff energy chosen somewhat arbitrarily, convergence checked in next program # # #
energy_cutoff = 300
energies = []
for n in n_k_points:
    if rank == 0:
        print 'Simulating ' + str(n) + 'kpoints...'
    mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)  # Recommended values for small systems
    calc = GPAW(mode=PW(energy_cutoff),  # use the LCAO basis mode
                h=0.18,  # grid spacing, recommended value in this course
                xc='PBE',  # Exchange-correlation functional
                mixer=mixer,  # Mixer
                # k-point grid - LOOP OVER LATER TO CHECK "CONVERGENCE"
                kpts=(n_k_points, n_k_points, n_k_points),
                txt='simulate_bulk_Al_GPAW.txt')  # name of GPAW output text file

    al_bulk.set_calculator(calc)
    total_energy = al_bulk.get_potential_energy()
    energies.append(total_energy)

with open('~/TIF035/HA2/bulk/1_converge_kpoints_bulk.txt', 'w') as textfile:
    texfile.write('number of k points, bulk_energy\n')
    for i in range(len(n_k_points)):
        textfile.write(str(n_k_points[i]) + ',' + str(energies[i]))
