from aluminium_nano_particle import aluminium_nano_particle

with open('~/TIF035/HA3/bulk/3_calc_optimal_lattice_spacing.txt', 'r') as textfile:
    lattice_parameter = float(next(textfile).split(':')[1])

N_x = 1
N_y = 1
N_z = 7

particle = aluminium_nano_particle()
particle.create_surface_Al('111', N_x, N_y, N_z, lattice_parameter)
n_k_points = [4, 5, 6, 7, 8, 9, 10]
energy_cutoff = 200
for n in n_k_points:
    particle.simulate_surface_Al('111', Nx, Ny, Nz, energy_cutoff, n)
    energies.append(particle.get_surface_energy('111'))

with open('~/TIF035/HA3/surface/4_converge_kpoints_bulk.txt', 'w') as textfile:
    texfile.write('number of k points, bulk_energy\n')
    for i in range(len(n_k_points)):
        textfile.write(str(n_k_points[i]) + ',' + str(energies[i]))
