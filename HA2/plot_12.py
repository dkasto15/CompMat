import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

''' Import files '''
filename1 = '1_converge_kpoints_bulk.txt'
filename2 = '2_converge_cutoff_energy_bulk.txt'

energy_1 = []
kpoints = []

energy_2 = []
energy_cutoff = []

with open('HA2/' + filename1, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        kpoints.append(float(line[0]))
        energy_1.append(float(line[1]))

with open('HA2/' + filename2, 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        energy_cutoff.append(float(line[0]))
        energy_2.append(float(line[1]))

''' Plotting '''
# if whichcase == 1:
#     x = kpoints
#     y = energy_1
# if whichcase == 2:
#     x = energy_cutoff
#     y = energy_2
fig = plt.figure()
ax_kpoints = fig.add_subplot(211)
ax_kpoints.minorticks_on()
ax_kpoints.grid(which='major', color='gray', linestyle='solid')
ax_kpoints.grid(which='minor', color='gray', linestyle='dashed')
ax_kpoints.plot(kpoints, energy_1)
ax_kpoints.set_ylabel('Simulated energy [atomic units]')
ax_kpoints.set_xlabel('Wave vector [atomic units]')
ax_kpoints.set_xlim([kpoints[0], kpoints[-1]])
ax_kpoints.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax_cutoff = fig.add_subplot(212)
ax_cutoff.minorticks_on()
ax_cutoff.grid(which='major', color='gray', linestyle='solid')
ax_cutoff.grid(which='minor', color='gray', linestyle='dashed')
ax_cutoff.plot(energy_cutoff, energy_2)
ax_cutoff.set_ylabel('Simulated energy [atomic units]')
ax_cutoff.set_xlabel('Cutoff energy [atomic units]')
ax_cutoff.set_xlim([energy_cutoff[0], energy_cutoff[-1]])
ax_cutoff.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

fig.subplots_adjust(hspace=0.5)

plt.show()
