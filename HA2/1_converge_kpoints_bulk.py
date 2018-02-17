from aluminium_nano_particle import aluminium_nano_particle

particle = aluminium_nano_particle()
particle.create_bulk_Al(particle.experimental_lattice_parameter)
n_k_points = [4, 5, 6, 7, 8, 9, 10]
energy_cutoff = 200
energies = []
for n in n_k_points:
    particle.simulate_bulk_Al(E_cutoff, n)
    energies.append(particle.get_bulk_energy())

with open('~/TIF035/HA3/bulk/1_converge_kpoints_bulk.txt', 'w') as textfile:
    texfile.write('number of k points, bulk_energy\n')
    for i in range(len(n_k_points)):
        textfile.write(str(n_k_points[i]) + ',' + str(energies[i]))
