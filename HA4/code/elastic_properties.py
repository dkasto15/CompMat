from eam_calculator import get_calc
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from ase.eos import EquationOfState
from ase.atoms import copy
from ase.build import bulk
from ase.phonons import Phonons
from ase.dft.kpoints import ibz_points, bandpath
global al_bulk


def main():
    # ''' Read in parameters for EAM calc '''
    # with open('HA4/results/fit_potential_output.txt', 'r') as textfile:
    #     line = next(textfile)
    #     line = line.split(',')
    #     A = float(line[0])
    #     lmbd = float(line[1])
    #     D = float(line[2])
    #     mu2 = float(line[3])

    # with open('HAlea')
    ''' Optimization parameters '''
    A = 1000  # eV
    lmbd = 3  # Å^(-1)
    D = 5  # Å
    mu2 = 1  # 2  # Å^(-1)
    param_0 = [A, lmbd, D, mu2]

    ''' Strains and stuff '''
    calc = get_calc((A, lmbd, D, mu2))
    e1 = 0.01
    e6 = 0.01
    energies = []

    # C11-C12
    al_bulk = bulk('Al', 'fcc', a=4.032)
    al_bulk.set_calculator(calc)
    ep_mat_1 = np.array([[e1, 0, 0], [0, -e1, 0], [0, 0, e1**2 / (1 - e1**2)]])
    energies.append(al_bulk.get_potential_energy())
    cell_0 = al_bulk.get_cell()
    al_bulk.set_cell(np.dot((np.eye(3) + ep_mat_1), cell_0))
    energies.append(al_bulk.get_potential_energy())
    V = al_bulk.get_volume()
    delta_E = energies[-1] - energies[0]
    C11_minus_C12 = delta_E / (V * e1**2)
    print(C11_minus_C12)

    # C11+C12
    al_bulk = bulk('Al', 'fcc', a=4.032)
    al_bulk.set_calculator(calc)
    e2 = e1
    ep_mat_12 = np.array([[e1, 0, 0], [0, e2, 0], [0, 0, 0]])
    energies.append(al_bulk.get_potential_energy())
    V = al_bulk.get_volume() # Equilibrium cell volume
    print(al_bulk.get_volume())
    cell_0 = al_bulk.get_cell()
    al_bulk.set_cell(np.dot((np.eye(3) + ep_mat_12), cell_0))
    print(al_bulk.get_volume())
    energies.append(al_bulk.get_potential_energy())
    delta_E = energies[-1] - energies[0]
    C11_plus_C12 = delta_E / (V * e1**2)
    print(C11_plus_C12)

    ## C11 and C12
    C11 = (C11_minus_C12 + C11_plus_C12) / 2
    C12 = C11_plus_C12 - C11

    # C44
    al_bulk = bulk('Al', 'fcc', a=4.032)
    al_bulk.set_calculator(calc)
    ep_mat_6 = np.array([[0, 0.5 * e6, 0], [0.5 * e6, 0, 0], [0, 0, e6**2 / (4 - e6**2)]])
    cell_0 = al_bulk.get_cell()
    al_bulk.set_cell(np.dot((np.eye(3) + ep_mat_6), cell_0))
    energies.append(al_bulk.get_potential_energy())
    V = al_bulk.get_volume()
    delta_E = energies[-1] - energies[0]
    C44 = 2 * delta_E / (V * e6**2)
    print(C44)

    B = (C11 + 2 * C12) / 3
    c_prim = (C11 * C12) / 2

    ''' Phonon calculator '''
    al_bulk = bulk('Al', 'fcc', a=4.032)
    N = 7
    ph = Phonons(al_bulk, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    # Read forces and assemble the dynamical matrix
    ph.read(acoustic=True, banana=True)

    # High-symmetry points in the Brillouin zone
    points = ibz_points['fcc']
    G = points['Gamma']
    X = points['X']
    W = points['W']
    K = points['K']
    L = points['L']
    U = points['U']

    point_names = ['$\Gamma$', 'X', 'U', 'L', '$\Gamma$', 'K']
    path = [G, X, U, L, G, K]

    # Band structure in meV
    path_kc, q, Q = bandpath(path, al_bulk.cell, 100)
    omega_kn = 1000 * ph.band_structure(path_kc)

    # Calculate phonon DOS
    omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=5e-4)
    omega_e *= 1000

    # Plot the band structure and DOS
    plt.figure(1, (8, 6))
    plt.axes([.1, .07, .67, .85])
    for n in range(len(omega_kn[0])):
        omega_n = omega_kn[:, n]
        plt.plot(q, omega_n, 'k-', lw=2)
    print('hej', np.sqrt(path_kc[:,0]**2 + path_kc[:,1]**2 + path_kc[:,2]**2 )/q)
    print(q)
    print(np.shape(q))
    print(np.shape(path_kc[:,0]))
    #
    # plt.xticks(Q, point_names, fontsize=18)
    # plt.yticks(fontsize=18)
    # plt.xlim(q[0], q[-1])
    # plt.ylabel("Frequency ($\mathrm{meV}$)", fontsize=22)
    # plt.grid('on')
    #
    # plt.axes([.8, .07, .17, .85])
    # plt.fill_between(dos_e, omega_e, y2=0, color='lightgrey', edgecolor='k', lw=1)
    # plt.ylim(0, 35)
    # plt.xticks([], [])
    # plt.yticks([], [])
    # plt.xlabel("DOS", fontsize=18)
    plt.show()

    ''' Sound velocity '''
    # point_names = ['$\Gamma$', 'X']
    path_100 = [G, X]

    # Band structure in meV
    # path_kc_100, q_100, Q_100 = bandpath(path_100, al_bulk.cell, 100) # Return list of k-points, list of x-coordinates and list of x-coordinates of special points.
    # print(path_kc_100)
    # print(q_100)
    # omega_kn_100 = 1000 * ph.band_structure(path_kc_100)
    # print('om', omega_kn_100)
    # print('om_test', test)
    #
    # plt.figure()
    # plt.plot(q_100, omega_kn_100)
    # plt.show()


if __name__ == '__main__':
    main()
