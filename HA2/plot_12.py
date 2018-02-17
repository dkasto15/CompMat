import matplot.pyplot as plt

energy = []
kpoints = []
with open('HA2/1_converge_kpoints_bulk.txt', 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        kpoints = float(line[0])
        energy = float(line[1])

fig = plt.figure()
ax_kpoints = fig.add_subplot(211)
ax_kpoints.plot(kpoints, energy)
ax_kpoints.set_ylabel('Simulated energy [atomic units]')
ax_kpoints.set_xlabel('Wave vector [atomic units]')
ax_kpoints.set_xlim([kpoints[0], kpoints[1]])
ax_kpoints.minor_ticks_on()
ax_kpoints.grid(True, which='both')

ax_energy = fig.add_subplot(212)

fig.show()
