from eam_calculator import get_calc
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from ase.eos import EquationOfState
from ase.atoms import copy
from ase.build import bulk
global n_sets, cells


def main():
    A = 1000  # eV
    lmbd = 3  # Å^(-1)
    D = 5  # Å
    mu2 = 2  # Å^(-1)
    x0 = [A, lmbd, D, mu2]

    file_tags = ['0.9', '1.0', '1.1']
    n_sets = len(file_tags)

    forces = []
    atoms_set = []
    ls_forces = True
    ls_energy = False
    ls_lattice = False

    atom_set = []  # Set of meassurements.
    for index, tag in enumerate(file_tags):
        atoms_set.append(read('HA4/downloads/snapshots_with_forces_xyz/res_POSCAR_' +
                              tag + '.xyz'))

        forces_tmp = atoms_set[index].get_forces()
        forces = np.hstack([forces, forces_tmp[:, 0], forces_tmp[:, 1], forces_tmp[:, 2]])

    forces = np.asarray(forces)

    a0 = 4.032
    E0 = -3.36
    al_bulk = bulk('Al', 'fcc', a=a0)

    ftol = 1e-8
    xtol = 1e-8
    gtol = 1e-8
    loss = 'linear'

    if ls_energy:
        res_cohesive_energy = least_squares(penalty_function_energies,
                                            x0,
                                            bounds=([0.5 * A, 0.5 * lmbd, 0.5 * D, 0.5 * mu2],
                                                    [2 * A, 2 * lmbd, 2 * D, 2 * mu2]),
                                            args=(al_bulk, E0),
                                            ftol=ftol,
                                            xtol=xtol,
                                            gtol=gtol,
                                            loss=loss,
                                            verbose=2)
        write_least_squares_output('energy', res_cohesive_energy)

    if ls_forces:
        res_forces = least_squares(penaly_function_forces,
                                   x0,
                                   bounds=([0.5 * A, 0.5 * lmbd, 0.5 * D, 0.5 * mu2],
                                           [1000000000, 2 * lmbd, 2 * D, 2 * mu2]),
                                   args=(atoms_set, forces),
                                   ftol=ftol,
                                   xtol=xtol,
                                   gtol=gtol,
                                   loss=loss,
                                   verbose=2)

        write_least_squares_output('forces', res_forces)

    if ls_lattice:
        res_lattice = least_squares(penalty_function_lattice_param,
                                    x0,
                                    args=(al_bulk, a0),
                                    ftol=ftol,
                                    xtol=xtol,
                                    gtol=gtol,
                                    loss=loss,
                                    verbose=2)

        write_least_squares_output('lattice_parameter', res_lattice)


def calc_forces_al(atoms_set, A, lmbd, D, mu2):
    forces = []
    calc = get_calc((A, lmbd, D, mu2))
    for atoms in atoms_set:
        atoms.set_calculator(calc)
        forces_tmp = atoms.get_forces()
        forces = np.hstack([forces, forces_tmp[:, 0], forces_tmp[:, 1], forces_tmp[:, 2]])
    return forces


def calc_cohesive_energy_al(al_bulk, A, lmbd, D, mu2):
    calc = get_calc((A, lmbd, D, mu2))
    al_bulk.set_calculator(calc)
    energy = al_bulk.get_potential_energy()
    return energy


def calc_lattice_parameter_al(al_bulk, A, lmbd, D, mu2):
    q = 0.03
    n_points = 14
    calc = get_calc((A, lmbd, D, mu2))
    # # # Find lattice constant with lowest energy # # #
    cell_0 = al_bulk.cell  # Unit cell object of the Al bulk
    al_bulk.set_calculator(calc)
    energies = []
    volumes = []

    for eps in np.linspace(-q, q, n_points):

        al_bulk.cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell

        # Calculate the potential energy for the Al bulk
        energies.append(al_bulk.get_potential_energy())
        volumes.append(al_bulk.get_volume())

    # Plot energies as a function of unit cell volume (directly related to latt. const.)
    eos = EquationOfState(volumes, energies)
    v0, E, B = eos.fit()
    # Latt. const. acc. to ASE doc., but why is this correct?
    a_calc = (4 * v0)**(1 / 3.0)
    al_bulk.cell = cell_0

    print('Lattice_parameter:', a_calc)
    return a_calc


def penaly_function_forces(x, atom_set, y):
    print('Forces: ' + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + ' ' + str(x[3]))
    return calc_forces_al(atom_set, x[0], x[1], x[2], x[3]) - y


def penalty_function_energies(x, al_bulk, y):
    print('Energy: ' + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + ' ' + str(x[3]))
    return calc_cohesive_energy_al(al_bulk, x[0], x[1], x[2], x[3]) - y


def penalty_function_lattice_param(x, al_bulk, y):
    print('Lattice: ' + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + ' ' + str(x[3]))
    return calc_lattice_parameter_al(al_bulk, x[0], x[1], x[2], x[3]) - y


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
