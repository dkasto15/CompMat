from eam_calculator import get_calc
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq
from ase.eos import EquationOfState
from ase.atoms import copy
from ase.build import bulk
global al_bulk


def main():
    ''' Location of input and experimental output data '''
    file_tags = ['0.9', '1.0', '1.1']

    ''' Optimization parameters '''
    A = 1000  # eV
    lmbd = 3  # Å^(-1)
    D = 5  # Å
    mu2 = 1  # 2  # Å^(-1)
    param_0 = [A, lmbd, D, mu2]

    ''' Input data '''
    al_bulk = bulk('Al', 'fcc', a=4.032)
    data_input = []  # Array for storing the different types of input data
    for index, tag in enumerate(file_tags):
        data_input.append(read('HA4/downloads/snapshots_with_forces_xyz/res_POSCAR_' +
                               tag + '.xyz'))
        print(data_input[index].get_calculator())
    data_input.append(al_bulk)

    ''' Observed output data '''
    a0 = 4.032
    E0 = -3.36
    data_exp = []   # Array for storing the experimental output data from the input data
    energy_evo = []  # Vector for storing the evolution of the energy vector
    lattice_evo = []  # Vector for storing the evolution of the lattice parameter

    for index, tag in enumerate(file_tags):
        forces = data_input[index].get_forces()
        data_exp = np.hstack([data_exp, forces[:, 0], forces[:, 1], forces[:, 2]])
    data_exp = np.append(data_exp, a0)
    data_exp = np.append(data_exp, E0)

    ftol = 1e-8
    xtol = 1e-8
    gtol = 1e-8
    loss = 'linear'

    w_force = 1
    w_E0 = 0  # 972
    w_a0 = 0  # 972
    weights = np.sqrt(np.array([w_force, w_E0, w_a0]))

    ''' Least squares optimization procedure '''
    res_cohesive_energy = leastsq(calc_residuals,
                                  param_0,
                                  # method='lm',
                                  args=(data_input, data_exp, weights),
                                  ftol=ftol,
                                  xtol=xtol,
                                  # diff_step=0.1,
                                  # gtol=gtol,
                                  # loss=loss,
                                  # verbose=2)
                                  )

    plt.plot(data_exp, data_sim)
    plt.show()


def calc_forces(data_input_forces, A, lmbd, D, mu2):
    calc = get_calc((A, lmbd, D, mu2))
    forces = []
    for atoms in data_input_forces:
        atoms.set_calculator(calc)
        forces_mat = atoms.get_forces()
        forces = np.hstack([forces, forces_mat[:, 0], forces_mat[:, 1], forces_mat[:, 2]])
    return forces


def calc_cohesive_energy(al_bulk_ASE, A, lmbd, D, mu2):
    calc = get_calc((A, lmbd, D, mu2))
    al_bulk_ASE.set_calculator(calc)
    return al_bulk_ASE.get_potential_energy()


def calc_lattice_parameter(al_bulk_ASE, A, lmbd, D, mu2):
    q = 0.2
    n_points = 15
    calc = get_calc((A, lmbd, D, mu2))

    # # # Find lattice constant with lowest energy # # #
    cell_0 = al_bulk_ASE.get_cell()  # Unit cell object of the Al bulk
    al_bulk_ASE.set_calculator(calc)
    energies = []
    volumes = []

    for eps in np.linspace(-q, q, n_points):

        al_bulk_ASE.set_cell((1 + eps) * cell_0)  # Adjust lattice constant of unit cell

        # Calculate the potential energy for the Al bulk
        energies.append(al_bulk_ASE.get_potential_energy())
        volumes.append(al_bulk_ASE.get_volume())

    #plt.plot(np.asarray(volumes)**(1 / 3), energies)
    # plt.show()
    # Plot energies as a function of unit cell volume (directly related to latt. const.)
    eos = EquationOfState(volumes, energies)
    v0, E, B = eos.fit()
    # Latt. const. acc. to ASE doc., but why is this correct?
    a_calc = v0**(1 / 3.0)
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
    data_sim = np.append(data_sim, calc_cohesive_energy(data_input[3], A, lmbd, D, mu2))

    print(data_sim[-3])
    print(data_sim[-2])
    print(data_sim[-1])

    w_force = weights[0]
    w_E0 = weights[1]
    w_a0 = weights[2]

    residuals = (data_sim - data_exp)
    residuals[:-2] = w_force * residuals[:-2]
    residuals[-2] = w_a0 * residuals[-2]
    residuals[-1] = w_E0 * residuals[-1]

    # print(residuals)
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

    else:
        f_least_squares.write(res.message)
    f_least_squares.close()


if __name__ == '__main__':
    main()
