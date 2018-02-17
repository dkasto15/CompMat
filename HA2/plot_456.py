import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
''' Import files '''
filename1 = '4_converge_kpoints_surface.txt'
filename2 = '5_converge_cutoff_energy_surface.txt'
filename3 = '6_converge_surface_depth.txt'

sigma_1 = []
kpoints = []

sigma_2 = []
energy_cutoff = []

sigma_3 = []
E_N = []
Nz = []

with open('HA2/' + filename1, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        kpoints.append(float(line[0]))
        sigma_1.append(float(line[2]))

with open('HA2/' + filename2, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        energy_cutoff.append(float(line[0]))
        sigma_2.append(float(line[2]))

with open('HA2/' + filename3, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        Nz.append(float(line[0]))
        E_N.append(float(line[0]))
        sigma_3.append(float(line[2]))

''' Plotting '''

fig = plt.figure()
ax_kpoints = fig.add_subplot(311)
ax_kpoints.minorticks_on()
ax_kpoints.grid(which='major', color='gray', linestyle='solid')
ax_kpoints.grid(which='minor', color='gray', linestyle='dashed')
ax_kpoints.plot(kpoints, sigma_1)
ax_kpoints.set_xlabel('Wave vector [atomic units]')
ax_kpoints.set_xlim([kpoints[0], kpoints[-1]])
ax_kpoints.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax_cutoff = fig.add_subplot(312)
ax_cutoff.minorticks_on()
ax_cutoff.grid(which='major', color='gray', linestyle='solid')
ax_cutoff.grid(which='minor', color='gray', linestyle='dashed')
ax_cutoff.plot(energy_cutoff, sigma_2)
ax_cutoff.set_ylabel('Simulated surface energy density [atomic units]')
ax_cutoff.set_xlabel('Cutoff energy [atomic units]')
ax_cutoff.set_xlim([energy_cutoff[0], energy_cutoff[-1]])
ax_cutoff.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax_Nz = fig.add_subplot(313)
ax_Nz.minorticks_on()
ax_Nz.grid(which='major', color='gray', linestyle='solid')
ax_Nz.grid(which='minor', color='gray', linestyle='dashed')
ax_Nz.plot(Nz, np.asarray(sigma_3), '.')
ax_Nz.set_xlabel('Atom depth in surface [atomic units]')
ax_Nz.set_xlim([Nz[0], Nz[-1]])
ax_Nz.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig.subplots_adjust(hspace=0.5)

plt.show()
