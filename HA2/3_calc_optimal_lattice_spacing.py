from aluminium_nano_particle import aluminium_nano_particle

particle = aluminium_nano_particle()
particle.create_bulk_aluminum(lattice_spacing=particle.a_al_experimental)

''' Calculates the lattice spacing that minimizes the
    potential energy for bulk aluminium. This is done by offsetting
    the experimentally derived lattice and then fitting the generated
    energies with GPAW's function EquationOfState.
'''
try:
    with open('~/TIF035/HA3/bulk/converged_kpoint.txt', 'r') as textfile:
        n_k_points = int(next(textfile))
    with open('~/TIF035/HA3/bulk/converged_cutoff_energy.txt', 'r') as textfile:
        E_cutoff = int(next(textfile))
except Exception as e:
    print('converged_kpoint.txt or converged_cutoff_energy.txt do not exist')
    exit(1)

q = 0.02  # Percentage difference from equalibrium for each point
n_points = 7  # n_points: Number of lattice constants to loop over to find equilibrium.
energies = []
volumes = []

# # # Find lattice constant with lowest energy # # #
cell_0 = particle.get_bulk_cell()  # Unit cell object of the Al bulk

for eps in np.linspace(-q, q, n_points):
    particle.set_cell((1 + eps) * cell_0)  # Adjust lattice constant of unit cell
    # Calculate the potential energy for the Al bulk
    simulate_bulk_Al(E_cutoff, n_k_points)
    energies.append(energy)
    volumes.append(volume)

    # Plot energies as a function of unit cell volume (directly related to latt. const.)
    eos = EquationOfState(volumes, energies)
    v0, E_bulk, B = eos.fit()
    eos.plot('Al_eos.png')

    # Latt. const. acc. to ASE doc., but why is this correct?
    a_calc = (4 * v0)**(1 / 3.0)

if rank == 0:
    print 'Energies: ' + str(energies)
    print 'Volumes: ' + str(volumes)
    print 'a_calc: ' + str(a_calc)

with open('~/TIF035/HA3/bulk/3_calc_optimal_lattice_spacing.txt', 'w') as textfile:
    textfile.write('Optimal lattice spacing: ' + str(a_calc))
    textfile.write('Corresponding energy: ' + str(E_bulk))
    textfile.write('Corresponding bulk modulus: ' + str(B))
    textfile.write('Energies, Volumes')
    for i in range(len(volumes)):
        textfile.write(str(energies[i]) + ',' + str(volumes[i]))
