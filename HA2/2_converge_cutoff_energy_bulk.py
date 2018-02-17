from aluminium_nano_particle import aluminium_nano_particle
try:
    with open('/bulk/converged_kpoint.txt', 'r') as textfile:
        n_k_points = int(next(textfile))
except Exception as e:
    print('converged_kpoint.txt does not exist')
    exit(1)

particle = aluminium_nano_particle()
particle.create_bulk_Al(particle.experimental_lattice_parameter)
cutoff_energies = [50, 100, 200, 300, 400, 500]
for cutoff_energy in cutoff_energies:
    particle.simulate_bulk_Al(cutoff_energy, n_k_points)
    energies.append(particle.get_bulk_energy())

with open('~/TIF035/HA3/bulk/2_converge_cutoff_energy_bulk.txt', 'w') as textfile:
    texfile.write('cutoff_energy, bulk_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' + str(energies[i]))
