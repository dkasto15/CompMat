import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy.misc import derivative


def main():
    V = []
    E_EAM = []

    with open('HA4/results/EV_FCC.txt', 'r') as textfile:
        next(textfile)
        for line in textfile:
            line = line.split(',')
            V.append(float(line[0]))
            E_EAM.append(float(line[1]))

    V = np.asarray(V)
    E_EAM = np.asarray(E_EAM)
    RMS_res = []

    n_vec = range(1, 10)

    ftol = 1e-5
    xtol = 1e-5
    gtol = 1e-5
    loss = 'linear'

    ''' Converge number of basis functions '''
    for n in n_vec:
        x_0 = np.ones(n)
        res = least_squares(residual, x_0, args=(V, E_EAM), verbose=0)
        RMS_res.append(np.sqrt(2 * res.cost / len(E_EAM)))

    fig_converge = plt.figure()
    ax_converge = fig_converge.add_subplot(111)
    ax_converge.plot(n_vec, RMS_res, color='blue')
    ax_converge.minorticks_on()
    ax_converge.grid(True, which='minor', linestyle='--')
    ax_converge.grid(True, which='major', linestyle='-')
    ax_converge.set_xlim(n_vec[0], n_vec[-1])
    ax_converge.set_xlabel('Number of basis functions [-]')
    ax_converge.set_ylabel('Root mean square of cost function [-]')

    ''' Using nominal number of basis functions '''
    n = 4
    x_0 = np.ones(n)
    res = least_squares(residual, x_0, args=(V, E_EAM), verbose=0)
    RMS_res.append(np.sqrt(2 * res.cost / len(E_EAM)))

    fig_compare = plt.figure()
    ax_compare = fig_compare.add_subplot(111)
    ax_compare.plot(V, E_EAM, '.', color='black')
    V = np.linspace(V[0], V[-1], 100)  # Adding more points for smoother plot
    ax_compare.plot(V, birch_energy(res.x, V), color='blue')
    ax_compare.set_xlim(V[0], V[-1])
    ax_compare.minorticks_on()
    ax_compare.grid(True, which='minor', linestyle='--')
    ax_compare.grid(True, which='major', linestyle='-')
    ax_compare.set_xlabel('Volume [Å$^{3}$]')
    ax_compare.set_ylabel('Cohesive energy [eV]')

    ''' Calculating bulk modulus using proposed equation '''
    atomic_units_to_newton = (82.387 / 51.421) * 1e-9
    angstrom_to_meter = 1e-10
    B = min(V) * birch_bulk(res.x, min(V))
    B_SI = B * atomic_units_to_newton / (angstrom_to_meter)**2
    B_GPa = B_SI / 1e9

    with open('HA4/results/Bulk_modulus_birch.txt', 'w') as textfile:
        textfile.write('Volume (Å^3): ' + str(min(V)) + '\n')
        textfile.write('Bulk modulus (GPa): ' + str(B_GPa) + '\n')
        textfile.write('Bulk modulus (a.u): ' + str(B))

    fig_converge.savefig("HA4/results/Converge_n_birch.png", bbox_inches='tight')
    fig_compare.savefig("HA4/results/Compare_data_birch.png", bbox_inches='tight')
    plt.show()


def residual(x, V, E_EAM):
    return birch_energy(x, V) - E_EAM


def birch_energy(x, V):
    E = 0
    for n in range(len(x)):
        E = E + x[n] * V**(-2 * n / 3)
    return E


def birch_bulk(x, V):
    B = 0
    dV = 1e-8
    ddE_ddV = (birch_energy(x, V + dV) - 2 * birch_energy(x, V) +
               birch_energy(x, V - dV)) / dV**2
    for n in range(len(x)):
        ddE_ddV = ((4 * n**2 + 6 * n) / 9 * V**((-2 * n - 6) / 3)) * x[n]
        B = B + V * ddE_ddV
    return B


if __name__ == '__main__':
    main()
