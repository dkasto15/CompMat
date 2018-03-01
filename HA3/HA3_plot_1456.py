import numpy as np
import matplotlib.pyplot as plt
wavefunction1 = []
wavefunction4 = []
wavefunction5 = []
wavefunction6 = []
r456 = []
r1 = []
with open('HA3/HA3_456.txt', 'r') as textfile:
    next(textfile)  # Skipping description
    next(textfile)  # Skipping names
    next(textfile)  # Skipping Energies
    next(textfile)  # Skipping eigenvalues

    for row in textfile:
        row = row.split(',')
        wavefunction4.append(float(row[0]))
        wavefunction5.append(float(row[1]))
        wavefunction6.append(float(row[2]))
        r456.append(float(row[3]))

wavefunction4 = np.asarray(wavefunction4)
wavefunction5 = np.asarray(wavefunction5)
wavefunction6 = np.asarray(wavefunction6)
r456 = np.asarray(r456)

with open('HA3/HA3_1.txt', 'r') as textfile:
    next(textfile)  # Skipping energy
    for row in textfile:
        row = row.split(',')
        wavefunction1.append(float(row[0]))
        r1.append(float(row[1]))

wavefunction1 = np.asarray(wavefunction1)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(r1, wavefunction1, '--', label='Hartree using basis with four s-orbitals')
ax.plot(r456, wavefunction4, '--', label='DFT with $V = V_H$')
ax.plot(r456, wavefunction5, '--', label='DFT with $V = V_H + V_x$')
ax.plot(r456, wavefunction6, '--', label='DFT with $V = V_H + V_x + V_c$')
ax.minorticks_on()
ax.grid(True, which='minor', linestyle='--')
ax.grid(True, which='major', linestyle='-')
ax.set_xlim([0, 5])
ax.set_ylim([0, 1.2 * max(wavefunction4)])
ax.set_xlabel('Radial coordinate [atomic units]')
ax.set_ylabel('Probability distribution function')
ax.set_title('Probability distribution function for Helium electrons achieved \n' +
             'four different computational methods')
ax.legend()

plt.show()
