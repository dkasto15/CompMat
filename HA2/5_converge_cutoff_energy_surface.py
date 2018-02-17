from aluminium_nano_particle import aluminium_nano_particle
try:
    with open('~/TIF035/HA3/surface/converged_kpoint.txt', 'r') as textfile:
        n_k_points = int(next(textfile))
    with open('~/TIF035/HA3/bulk/3_calc_optimal_lattice_spacing.txt', 'r') as textfile:
        lattice_parameter = float(next(textfile).split(':')[1])
except Exception as e:
    print('converged_kpoint.txt does not exist')
    exit(1)


N_x = 1
N_y = 1
N_z = 7

particle = aluminium_nano_particle()
particle.create_surface_Al('111', N_x, N_y, N_z, lattice_parameter)
cutoff_energies = [50, 100, 200, 300, 400, 500]
for cutoff_energy in cutoff_energies:
    particle.simulate_surface_Al('111', Nx, Ny, Nz, cutoff_energy, n_k_points)
    energies.append(particle.get_surface_energy('111'))

particle.create_surface_Al('100', N_x, N_y, N_z, lattice_parameter)
particle.simulate_surface_Al(cutoff_energies[-1], n_k_points)
energy100 = particle.get_surface_energy('100')

with open('~/TIF035/HA3/surface/5_converge_cutoff_energy_surface.txt', 'w') as textfile:
    texfile.write('cutoff_energy, bulk_energy\n')
    for i in range(len(cutoff_energies)):
        textfile.write(str(cutoff_energies[i]) + ',' + str(energies[i]))
