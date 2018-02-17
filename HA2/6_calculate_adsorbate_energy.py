from aluminium_nano_particle import aluminium_nano_particle

try:
    with open('~/TIF035/HA3/surface/converged_kpoint.txt', 'r') as textfile:
        n_k_points = int(next(textfile))
    with open('~/TIF035/HA3/surface/converged_cutoff_energy.txt', 'r') as textfile:
        cutoff_energy = int(next(textfile))
    with open('~/TIF035/HA3/bulk/3_calc_optimal_lattice_spacing.txt', 'r') as textfile:
        lattice_parameter = float(next(textfile).split(':')[1])
    with open('~/TIF035/HA3/bulk/bulk_energy.txt', 'r') as textfile:
        energy_bulk = float(next(textfile).split(':')[1])
except Exception as e:
    print('converged_kpoint.txt does not exist')
    exit(1)


N_x = 1
N_y = 1
N_z = 7

particle = aluminium_nano_particle()

particle.create_surface_Al('111', N_x, N_y, N_z, lattice_parameter)
particle.create_surface_Al('100', N_x, N_y, N_z, lattice_parameter)

energy_111 = particle.simulate_surface_Al('111', Nx, Ny, Nz, cutoff_energy, n_k_points)
energy_100 = particle.simulate_surface_Al('100', Nx, Ny, Nz, cutoff_energy, n_k_points)

particle.add_CO_adsorbate('111', 5)
particle.add_CO_adsorbate('100', 5)

energy_CO = particle.simulate_CO(cutoff_energy, n_k_points)

particle.fixate_surface_atoms('111')
particle.fixate_surface_atoms('100')

particle.relax_CO_adsorbate('111')
particle.relax_CO_adsorbate('100')

energy_111_CO = particle.simulate_surface_Al('111', Nx, Ny, Nz, cutoff_energy, n_k_points)
energy_100_CO = particle.simulate_surface_Al('100', Nx, Ny, Nz, cutoff_energy, n_k_points)

area_111 = particle.get_surface_area('111')
sigma_111 = (1 / (2.0 * area111)) * (energy_111 - N_x * N_y * energy_bulk)

area_100 = particle.get_surface_area('100')
sigma_100 = (1 / (2.0 * area100)) * (energy_100 - N_x * N_y * energy_bulk)

sigma_111_CO = sigma_111 + theta * (energy_111_CO - energy_111 - energy_CO) / area_111
sigma_100_CO = sigma_100 + theta * (energy_100_CO - energy_100 - energy_CO) / area_100


with open('~/TIF035/HA3/surface/surface_energies_Al.txt', 'w') as textfile:
    textfile.write('111:' + str(energies))
    textfile.write('100:' + str(energy100))
