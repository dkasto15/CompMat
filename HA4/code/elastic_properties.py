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
from mpl_toolkits.mplot3d import Axes3D
from ase.calculators.emt import EMT
global al_bulk


def main():
    ''' Read in parameters for EAM calc '''
    with open('HA4/results/fit_potential_output_full.txt', 'r') as textfile:
        line = next(textfile)
        line = line.split(',')
        A = float(line[0])
        lmbd = float(line[1])
        D = float(line[2])
        mu2 = float(line[3])

    # with open('HAlea')
    # ''' Optimization parameters '''
    # A = 1000  # eV
    # lmbd = 3  # Å^(-1)
    # D = 5  # Å
    # mu2 = 1  # 2  # Å^(-1)
    # param_0 = [A, lmbd, D, mu2]

    ''' Strains and stuff '''
    calc = get_calc((A, lmbd, D, mu2))
    # calc = EMT()
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
    V = 4 * al_bulk.get_volume()
    delta_E = energies[-1] - energies[0]
    C11_minus_C12 = delta_E / (V * e1**2)
    print(C11_minus_C12)

    # C11+C12
    al_bulk = bulk('Al', 'fcc', a=4.032)
    al_bulk.set_calculator(calc)
    e2 = e1
    ep_mat_12 = np.array([[e1, 0, 0], [0, e2, 0], [0, 0, 0]])
    energies.append(al_bulk.get_potential_energy())
    V = 4 * al_bulk.get_volume() # Equilibrium cell volume
    cell_0 = al_bulk.get_cell()
    al_bulk.set_cell(np.dot((np.eye(3) + ep_mat_12), cell_0))
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
    V = 4 * al_bulk.get_volume()
    delta_E = energies[-1] - energies[0]
    C44 = 2 * delta_E / (V * e6**2)
    print(C44)

    B = (C11 + 2 * C12) / 3
    print('B', B)
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

    # # Check band path
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot(path_kc[:,0], path_kc[:,1], path_kc[:,2])
    # plt.show()

    # Calculate phonon DOS
    # omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=5e-4)
    # omega_e *= 1000
    #
    # # Plot the band structure and DOS
    # plt.figure(1, (8, 6))
    # plt.axes([.1, .07, .67, .85])
    # for n in range(len(omega_kn[0])):
    #     omega_n = omega_kn[:, n]
    #     plt.plot(q, omega_n, 'k-', lw=2)
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
    # plt.show()

    ''' Sound velocity '''
    # point_names = ['$\Gamma$', 'X']
    path_100 = [G, X]

    # Band structure in meV
    path_kc_100, q_100, Q_100 = bandpath(path_100, al_bulk.cell, 100) # Return list of k-points, list of x-coordinates and list of x-coordinates of special points.
    omega_kn_100 = 1000 * ph.band_structure(path_kc_100)

    # # Find the longitudinal curve (the one that is not initially overlapping)
    # print(omega_kn_100[0:10,0])
    # print(omega_kn_100[0:10,1])
    # print(omega_kn_100[0:10,2]) # <-- This one!

    k = np.sqrt(path_kc_100[:,0]**2 + path_kc_100[:,1]**2 + path_kc_100[:,2]**2) # [Å^-1]
    convert_meV_to_1_over_s = (1 / 1000) * (1 / (6.582119514 * 10**(-16)))
    # print(omega_kn_100[1,2])

    omega_long_at_q_to_0 = omega_kn_100[1,2] * convert_meV_to_1_over_s
    # omega_long_at_q_to_1 = omega_kn_100[2,2] * convert_meV_to_1_over_s

    c_s = omega_long_at_q_to_0 * 10**(-10) / k[1] # Speed of sound, [m/s]
    # c_s = 10**(-10) * ((omega_long_at_q_to_1 - omega_long_at_q_to_0)  / (k[2] - k[1])) # Speed of sound, [m/s]
    print(c_s)
    #
    # convert_u_to_kg = 1.66054 * 10**(-27)
    # convert_kg_to_eV_c2 = (2.99792 * 10**8)**2
    # m_Al = al_bulk.get_masses()[0] * convert_u_to_kg * convert_kg_to_eV_c2 # [eV * (s^2/m^2)]
    # nbr_of_atoms_UC = 4 # Number of atoms per unit cell for fcc
    # V_Al = nbr_of_atoms_UC * al_bulk.get_volume() # [Å^3]
    # rho_Al = m_Al * nbr_of_atoms_UC / V_Al # [eV * (s^2/m^2) / Å^3]
    # young = c_s**2 * rho_Al
    #
    # print(C11)
    # print(young)


    plt.figure()
    plt.plot(q_100, omega_kn_100)
    plt.show()


if __name__ == '__main__':
    main()
