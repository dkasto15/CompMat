import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
''' Import files '''
filename1 = '8_converge_cutoff_CO.txt'

energy_1 = []
energy_cutoff = []

with open('HA2/' + filename1, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        energy_cutoff.append(float(line[0]))
        energy_1.append(float(line[1]))

''' Plotting '''
# if whichcase == 1:
#     x = kpoints
#     y = energy_1
# if whichcase == 2:
#     x = energy_cutoff
#     y = energy_2
fig = plt.figure()
ax_kpoints = fig.add_subplot(111)
ax_kpoints.minorticks_on()
ax_kpoints.grid(which='major', color='gray', linestyle='solid')
ax_kpoints.grid(which='minor', color='gray', linestyle='dashed')
ax_kpoints.plot(energy_cutoff, energy_1)
ax_kpoints.set_ylabel('Simulated energy [eV]')
ax_kpoints.set_xlabel('Cutoff energy [eV]')
ax_kpoints.set_xlim([energy_cutoff[0], energy_cutoff[-1]])
ax_kpoints.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

plt.show()
fig.savefig('converge_CO.eps', bbox_inches='tight')
