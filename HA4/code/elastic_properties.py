from eam_calculator import get_calc
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from ase.eos import EquationOfState
from ase.atoms import copy
from ase.build import bulk
from ase.phonons import Phonons
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
    energies = []
    volumes = []
    al_bulk_ASE.set_calculator(calc)

    # # # Find lattice constant with lowest energy # # #
    cell_0 = al_bulk_ASE.get_cell()  # Unit cell object of the Al bulk

    # DIRECTION
    cell_orig = cell_0

    E1 = al_bulk_ASE.get_potential_energy()
    al_bulk_ASE.set_cell((1 + eps) * cell_0)
    E2 = al_bulk_ASE.get_potential_energy()

    direction = (E2 - E1) / abs((E2 - E1))
    print('Direction', direction)

    # Go one step in the right
    al_bulk_ASE.set_cell(cell_orig)
    energies.append(al_bulk_ASE.get_potential_energy())
    volumes.append(al_bulk_ASE.get_volume())

    cell_0 = al_bulk_ASE.get_cell()

    al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)  # Adjust lattice constant of unit cell

    energies.append(al_bulk_ASE.get_potential_energy())
    volumes.append(al_bulk_ASE.get_volume())

    # Go two steps back
    while (energies[-1] - energies[-2]) < 0:
        cell_0 = al_bulk_ASE.get_cell()
        # Adjust lattice constant of unit cell
        al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)

        # Calculate the potential energy for the Al bulk
        energies.append(al_bulk_ASE.get_potential_energy())
        volumes.append(al_bulk_ASE.get_volume())

    cell_0 = al_bulk_ASE.get_cell()
    al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)  # Adjust lattice constant of unit cell

    # Calculate the potential energy for the Al bulk
    energies.append(al_bulk_ASE.get_potential_energy())
    volumes.append(al_bulk_ASE.get_volume())
    print(energies)
    print(volumes)
    #plt.plot((4 * np.asarray(volumes))**(1 / 3), energies)
    # plt.show()
    # Plot energies as a function of unit cell volume (directly related to latt. const.)
    eos = EquationOfState(volumes, energies)
    v0, E, B = eos.fit()
    # Latt. const. acc. to ASE doc., but why is this correct?
    a_calc = (4 * v0)**(1 / 3.0)
    # a_calc = v0**(1 / 3.0)
    print('a=' + str(a_calc))

    al_bulk_ASE.set_cell(cell_0)

    return a_calc


def calc_residuals(optimization_params, data_input, data_exp, weights):
    print('Optimization params: ', optimization_params)
    A = optimization_params[0]
    lmbd = optimization_params[1]
    D = optimization_params[2]
    mu2 = optimization_params[3]

    data_input_forces = data_input[0:3]
    data_sim = calc_forces(data_input_forces, A, lmbd, D, mu2)
    data_sim = np.append(data_sim, calc_lattice_parameter(data_input[3], A, lmbd, D, mu2))
    data_sim = np.append(data_sim, calc_cohesive_energy(data_input[4], A, lmbd, D, mu2))

    w_force = weights[0]
    w_E0 = weights[1]
    w_a0 = weights[2]

    residuals = (data_sim - data_exp)
    residuals[:-2] = w_force * residuals[:-2]
    residuals[-2] = w_a0 * residuals[-2]
    residuals[-1] = w_E0 * residuals[-1]

    return residuals


def write_least_squares_output(name, res):
    f_least_squares = open('HA4/results/least_squares_' + name + '.txt', 'w')
    if res.success:
        f_least_squares.write('Success!\n')
        f_least_squares.write(res.message + '\n')

        text = 'Output parameters: '
        for param in res.x:
            text = text + str(param) + ' '
        text = text[:-1]
        text = text + '\n'

        text = text + 'Cost: ' + str(res.cost) + '\n'

        text = text + 'Residuals: '
        for residual in res.fun:
            text = text + str(residual) + ' '
        text = text[:-1]
        text = text + '\n'

        text = text + 'Gradient of cost function at solution: ' + str(res.grad) + '\n'
        text = text + 'Optimality at solution: ' + str(res.optimality) + '\n'
        text = text + 'Number of function evaluations: ' + str(res.nfev) + '\n'
        text = text + 'Number of jacobian evaluations: ' + str(res.njev)
        f_least_squares.write(text)

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
    print(al_bulk.get_volume())
    cell_0 = al_bulk.get_cell()
    al_bulk.set_cell(np.dot((np.eye(3) + ep_mat_12), cell_0))
    print(al_bulk.get_volume())
    energies.append(al_bulk.get_potential_energy())
    V = al_bulk.get_volume()
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


#     ''' Observed output data '''
#     a0 = 4.032
#     E0 = -3.36
#
#     with open('HA4/results/output.txt', 'a') as textfile:
#         line = ''
#         for el in res.x:
#             line = line + str(el) + ','
#         line = line[:-1]
#         textfile.write('\n')
#         line = line + 'Lattice parameter: ' + str(a0_final) + '\n'
#         line = line + 'Cohesive energy: ' + str(E0_final) + '\n\n'
#
#
# def calc_forces(data_input_forces, A, lmbd, D, mu2):
#     calc = get_calc((A, lmbd, D, mu2))
#     forces = []
#     for atoms in data_input_forces:
#         atoms.set_calculator(calc)
#         forces_mat = atoms.get_forces()
#         forces = np.hstack([forces, forces_mat[:, 0], forces_mat[:, 1], forces_mat[:, 2]])
#     return forces
#
#
# def calc_cohesive_energy(al_bulk_ASE, A, lmbd, D, mu2):
#     calc = get_calc((A, lmbd, D, mu2))
#     al_bulk_ASE.set_calculator(calc)
#     energy = al_bulk_ASE.get_potential_energy()
#     print('Cohesive energy: ', energy)
#     return energy
#
#
# def calc_lattice_parameter(al_bulk_ASE, A, lmbd, D, mu2):
#     calc = get_calc((A, lmbd, D, mu2))
#     eps = 0.001
#     energies = []
#     volumes = []
#     al_bulk_ASE.set_calculator(calc)
#
#     # # # Find lattice constant with lowest energy # # #
#     cell_0 = al_bulk_ASE.get_cell()  # Unit cell object of the Al bulk
#
#     # DIRECTION
#     cell_orig = cell_0
#
#     E1 = al_bulk_ASE.get_potential_energy()
#     al_bulk_ASE.set_cell((1 + eps) * cell_0)
#     E2 = al_bulk_ASE.get_potential_energy()
#
#     direction = (E2 - E1) / abs((E2 - E1))
#     print('Direction', direction)
#
#     # Go one step in the right
#     al_bulk_ASE.set_cell(cell_orig)
#     energies.append(al_bulk_ASE.get_potential_energy())
#     volumes.append(al_bulk_ASE.get_volume())
#
#     cell_0 = al_bulk_ASE.get_cell()
#
#     al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)  # Adjust lattice constant of unit cell
#
#     energies.append(al_bulk_ASE.get_potential_energy())
#     volumes.append(al_bulk_ASE.get_volume())
#
#     # Go two steps back
#     while (energies[-1] - energies[-2]) < 0:
#         cell_0 = al_bulk_ASE.get_cell()
#         # Adjust lattice constant of unit cell
#         al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)
#
#         # Calculate the potential energy for the Al bulk
#         energies.append(al_bulk_ASE.get_potential_energy())
#         volumes.append(al_bulk_ASE.get_volume())
#
#     cell_0 = al_bulk_ASE.get_cell()
#     al_bulk_ASE.set_cell((1 - direction * eps) * cell_0)  # Adjust lattice constant of unit cell
#
#     # Calculate the potential energy for the Al bulk
#     energies.append(al_bulk_ASE.get_potential_energy())
#     volumes.append(al_bulk_ASE.get_volume())
#     print(energies)
#     print(volumes)
#     #plt.plot((4 * np.asarray(volumes))**(1 / 3), energies)
#     # plt.show()
#     # Plot energies as a function of unit cell volume (directly related to latt. const.)
#     eos = EquationOfState(volumes, energies)
#     v0, E, B = eos.fit()
#     # Latt. const. acc. to ASE doc., but why is this correct?
#     a_calc = (4 * v0)**(1 / 3.0)
#     # a_calc = v0**(1 / 3.0)
#     print('a=' + str(a_calc))
#
#     al_bulk_ASE.set_cell(cell_0)
#
#     return a_calc
#
#
# def calc_residuals(optimization_params, data_input, data_exp, weights):
#     print('Optimization params: ', optimization_params)
#     A = optimization_params[0]
#     lmbd = optimization_params[1]
#     D = optimization_params[2]
#     mu2 = optimization_params[3]
#
#     data_input_forces = data_input[0:3]
#     data_sim = calc_forces(data_input_forces, A, lmbd, D, mu2)
#     data_sim = np.append(data_sim, calc_lattice_parameter(data_input[3], A, lmbd, D, mu2))
#     data_sim = np.append(data_sim, calc_cohesive_energy(data_input[4], A, lmbd, D, mu2))
#
#     # print(data_sim[-3])
#     # print(data_sim[-2])
#     # print(data_sim[-1])
#
#     w_force = weights[0]
#     w_E0 = weights[1]
#     w_a0 = weights[2]
#
#     residuals = (data_sim - data_exp)
#     residuals[:-2] = w_force * residuals[:-2]
#     residuals[-2] = w_a0 * residuals[-2]
#     residuals[-1] = w_E0 * residuals[-1]
#
#     # print(residuals)
#     return residuals
#
#
# def write_least_squares_output(name, res):
#     f_least_squares = open('HA4/results/least_squares_' + name + '.txt', 'w')
#     if res.success:
#         f_least_squares.write('Success!\n')
#         f_least_squares.write(res.message + '\n')
#
#         text = 'Output parameters: '
#         for param in res.x:
#             text = text + str(param) + ' '
#         text = text[:-1]
#         text = text + '\n'
#
#         text = text + 'Cost: ' + str(res.cost) + '\n'
#
#         text = text + 'Residuals: '
#         for residual in res.fun:
#             text = text + str(residual) + ' '
#         text = text[:-1]
#         text = text + '\n'
#
#         text = text + 'Gradient of cost function at solution: ' + str(res.grad) + '\n'
#         text = text + 'Optimality at solution: ' + str(res.optimality) + '\n'
#         text = text + 'Number of function evaluations: ' + str(res.nfev) + '\n'
#         text = text + 'Number of jacobian evaluations: ' + str(res.njev)
#         f_least_squares.write(text)
#
#     else:
#         f_least_squares.write(res.message)
#     f_least_squares.close()


if __name__ == '__main__':
    main()
