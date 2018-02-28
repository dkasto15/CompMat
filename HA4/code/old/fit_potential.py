from eam_calculator import get_calc
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from ase.eos import EquationOfState
from ase.atoms import copy

global n_sets, cells

A = 100  # eV
lmbd = 3  # Å^(-1)
D = 5  # Å
mu2 = 2  # Å^(-1)
x0 = [A, lmbd, D, mu2]

file_tags = ['0.9', '1.0', '1.1']
n_sets = len(file_tags)

forces = []
cohesive_energies = []
lattice_parameters = []
ls_forces = True
ls_energy = True
ls_lattice = True
cells = []
atom_set = []  # Set of meassurements.
for index, tag in enumerate(file_tags):
    atom_set.append(read('HA4/downloads/snapshots_with_forces_xyz/res_POSCAR_' +
                         tag + '.xyz'))
    forces_tmp = atom_set[index].get_forces()
    cohesive_energy_tmp = atom_set[index].get_potential_energy()
    lattice_param_tmp = (atom_set[index].get_volume())**(1 / 3)  # Lattice spaing 13!?!?
    cells_tmp = atom_set[index].cell

    cells.append(cells_tmp)
    lattice_parameters.append(lattice_param_tmp)
    forces = np.hstack([forces, forces_tmp[:, 0], forces_tmp[:, 1], forces_tmp[:, 2]])
    cohesive_energies.append(cohesive_energy_tmp)

forces = np.asarray(forces)
cohesive_energies = np.asarray(cohesive_energies)
lattice_parameters = np.asarray(lattice_parameters)


a0 = lattice_parameters
ftol = 1e-8
xtol = 1e-8
gtol = 1e-8
loss = 'linear'


def calc_forces_al(atoms_set, A, lmbd, D, mu2):
    forces = []
    calc = get_calc((A, lmbd, D, mu2))
    for atoms in atom_set:
        atoms.set_calculator(calc)
        forces_tmp = atoms.get_forces()
        forces = np.hstack([forces, forces_tmp[:, 0], forces_tmp[:, 1], forces_tmp[:, 2]])
    return forces


def calc_cohesive_energy_al(atom_set, A, lmbd, D, mu2):
    global n_sets
    energy_set = np.zeros(n_sets)
    calc = get_calc((A, lmbd, D, mu2))
    for index, atoms in enumerate(atom_set):
        atoms.set_calculator(calc)
        energy_set[index] = atoms.get_potential_energy()
    return energy_set


def calc_lattice_parameter_al(atom_set, A, lmbd, D, mu2):
    global n_sets, cells
    q = 0.2
    n_points = 10
    lattice_param_set = np.zeros(n_sets)
    calc = get_calc((A, lmbd, D, mu2))

    for index, atoms in enumerate(atom_set):
        # # # Find lattice constant with lowest energy # # #
        cell_0 = atoms.cells[index]  # Unit cell object of the Al bulk
        atoms.set_calculator(calc)
        energies = []
        volumes = []
        inval = np.linspace(-q, q, n_points)
        for eps in inval:

            atoms.cell = (1 + eps) * cell_0  # Adjust lattice constant of unit cell

            # Calculate the potential energy for the Al bulk
            energies.append(atoms.get_potential_energy())
            volumes.append(atoms.get_volume())

        # Plot energies as a function of unit cell volume (directly related to latt. const.)
        eos = EquationOfState(volumes, energies)
        v0, E, B = eos.fit()
        # Latt. const. acc. to ASE doc., but why is this correct?
        a_calc = v0**(1 / 3.0)

        n_min = np.argmin(energies)
        atoms.cell = (1 + inval[n_min]) * cell_0
        lattice_param_set[index] = a_calc
    return lattice_param_set


def penaly_function_forces(x, atom_set, y):
    print('Forces: ' + str(x[0]) + ' ' + str(x[1]) + ' ' str(x[2]) + ' ' str(x[3]))
    return calc_forces_al(atom_set, x[0], x[1], x[2], x[3]) - y


def penalty_function_energies(x, atom_set, y):
    print('Energy: ' + str(x[0]) + ' ' + str(x[1]) + ' ' str(x[2]) + ' ' str(x[3]))
    return calc_cohesive_energy_al(atom_set, x[0], x[1], x[2], x[3]) - y


def penalty_function_lattice_param(x, atom_set, y):
    print('Lattice: ' + str(x[0]) + ' ' + str(x[1]) + ' ' str(x[2]) + ' ' str(x[3]))
    return calc_lattice_parameter_al(atom_set, x[0], x[1], x[2], x[3]) - y


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


if ls_energy:
    res_cohesive_energy = least_squares(penalty_function_energies,
                                        x0,
                                        args=(atom_set, cohesive_energies),
                                        ftol=ftol,
                                        xtol=xtol,
                                        gtol=gtol,
                                        loss=loss,
                                        verbose=2)
    write_least_squares_output('energy', res_cohesive_energy)

if ls_forces:
    res_forces = least_squares(penaly_function_forces,
                               x0,
                               args=(atom_set, forces),
                               ftol=ftol,
                               xtol=xtol,
                               gtol=gtol,
                               loss=loss,
                               verbose=2)

    write_least_squares_output('forces', res_forces)

if ls_lattice:
    res_lattice = least_squares(penalty_function_lattice_param,
                                x0,
                                args=(atom_set, lattice_parameters),
                                ftol=ftol,
                                xtol=xtol,
                                gtol=gtol,
                                loss=loss,
                                verbose=2)

    write_least_squares_output('lattice_parameter', res_lattice, f_least_squares)

with open('HA4/results/least_squares_solutions.txt', 'w') as textfile:
    textfile.write('cohesive energy, force, lattice parameter')
    for i in range(len(x0)):
        a = res_cohesive_energy.x[i]
        b = res_forces.x[i]
        c = res_lattice.x[i]
        textfile.write(str(a) + ', ' + str(b), ', ' + str(c))
