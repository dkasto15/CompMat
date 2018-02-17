import matplotlib.pyplot as plt

energy_kpoints = []
kpoints = []
with open('HA2/1_converge_kpoints_bulk.txt', 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        kpoints.append(float(line[0]))
        energy.append(float(line[1]))

fig = plt.figure()
ax_kpoints = fig.add_subplot(211)
ax_kpoints.minorticks_on()
ax_kpoints.grid(which='major', color='gray', linestyle='solid')
ax_kpoints.grid(which='minor', color='gray', linestyle='dashed')
ax_kpoints.plot(kpoints, energy)
ax_kpoints.set_ylabel('Simulated energy [atomic units]')
ax_kpoints.set_xlabel('Wave vector [atomic units]')
ax_kpoints.set_xlim([kpoints[0], kpoints[-1]])

ax_energy = fig.add_subplot(212)
ax_cutoff.minorticks_on()
ax_cutoff.grid(which='major', color='gray', linestyle='solid')
ax_cutoff.grid(which='minor', color='gray', linestyle='dashed')
ax_cutoff.plot(cutoff, energy)
ax_cutoff.set_ylabel('Simulated energy [atomic units]')
ax_cutoff.set_xlabel('Wave vector [atomic units]')
ax_cutoff.set_xlim([cutoff[0], cutoff[-1]])

fig.subplots_adjust(hspace=0.5)

plt.show()
