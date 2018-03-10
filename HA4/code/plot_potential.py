import matplotlib.pyplot as plt
import numpy as np


def pair_potential(r, A, lmbd):
    return A * np.exp(-lmbd * r)


def embedding_function(rho, D):
    return -D * np.sqrt(rho)


def electron_density(r, mu2):
    return np.exp(-mu2 * r)


def cutoff_function(x):
    fc = np.zeros(np.shape(x))
    for i in range(len(x)):
        if x[i] <= 0:
            fc[i] = x[i]**4 / (1 + x[i]**4)
        else:
            fc[i] = 0
    return fc


with open('HA4/results/fit_potential_output.txt', 'r') as textfile:
    line = next(textfile)
    line = line.split(',')
    A = float(line[0])
    lmbd = float(line[1])
    D = float(line[2])
    mu2 = float(line[3])
    next(textfile)
    next(textfile)
    line = next(textfile)
    line = line.split(':')
    RMS_res = float(line[1])

forces_eam = []
forces_dft = []

with open('HA4/results/fit_potential_force_components.txt', 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        forces_eam.append(float(line[0]))
        forces_dft.append(float(line[1]))
forces_eam = np.asarray(forces_eam)
forces_dft = np.asarray(forces_dft)

r_c = 6.5
r = np.linspace(0, r_c + 0.5, 100)

rho = electron_density(r, mu2) * cutoff_function((r - r_c) / 3)
F = embedding_function(rho, D)
V = pair_potential(r, A, lmbd)


''' Plotting '''
y_axis_coords = (-0.1, 0.5)
fig_EAM_potential = plt.figure(figsize=(6, 8))
ax_density = fig_EAM_potential.add_subplot(311)
ax_embedding = fig_EAM_potential.add_subplot(312)
ax_pair_potential = fig_EAM_potential.add_subplot(313)

ax_density.plot(r, rho, color='blue')
ax_density.set_xlim(0, r[-1])
ax_density.minorticks_on()
ax_density.grid(True, which='minor', linestyle='--')
ax_density.grid(True, which='major', linestyle='-')
ax_density.set_ylabel('Electron density [Å$^{-3}$]')
ax_density.yaxis.set_label_coords(y_axis_coords[0], y_axis_coords[1])

ax_embedding.plot(r, F, color='blue')
ax_embedding.set_xlim(0, r[-1])
ax_embedding.minorticks_on()
ax_embedding.grid(True, which='minor', linestyle='--')
ax_embedding.grid(True, which='major', linestyle='-')
ax_embedding.set_ylabel('Embedding function [eV]')
ax_embedding.yaxis.set_label_coords(y_axis_coords[0], y_axis_coords[1])

ax_pair_potential.plot(r, V, color='blue')
ax_pair_potential.set_xlim(0, r[-1])
ax_pair_potential.minorticks_on()
ax_pair_potential.grid(True, which='minor', linestyle='--')
ax_pair_potential.grid(True, which='major', linestyle='-')
ax_pair_potential.set_ylabel('Pair potential [eV]')
ax_pair_potential.set_xlabel('Radial coordinate [Å]')
ax_pair_potential.yaxis.set_label_coords(y_axis_coords[0], y_axis_coords[1])

fig_compare_forces = plt.figure()
ax_compare_forces = fig_compare_forces.add_subplot(111)
ax_compare_forces.plot(forces_dft, forces_eam, '.', color='blue',
                       label='RMS of residual: ' + str(np.round(RMS_res * 1000) / 1000))
ax_compare_forces.minorticks_on()
ax_compare_forces.grid(True, which='minor', linestyle='--')
ax_compare_forces.grid(True, which='major', linestyle='-')
ax_compare_forces.set_ylabel('Forces EAM [atomic units]')
ax_compare_forces.set_xlabel('Forces DFT [atomic units]')
ax_compare_forces.legend()

fig_EAM_potential.savefig("HA4/results/EAM_potential.png", bbox_inches='tight')
fig_compare_forces.savefig("HA4/results/Compare_forces.png", bbox_inches='tight')


plt.show()
