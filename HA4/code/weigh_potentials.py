x_energy = []
x_force = []
x_lattice = []


with open('HA4/results/least_squares_solutions.txt', 'r') as textfile:
    next(textfile)
    for line in textfile:
        line = line.split(',')
        x_energy.append(float(line[0]))
        x_force.append(float(line[1]))
        x_lattice.append(float(line[2]))

global w_force
w_force = 1


def chi_squared(w_energy, w_force, w_lattice, cost_energy, cost_force, cost_lattice):
    return w_force * cost_force + w_energy * cost_energy + w_lattice * cost_lattice


def min_fun(x, cost, y):
    global w_force
    return chi_squared(w_force, x[0], x[1], cost[0], cost[1], cost[2])


'''TODO: APPLY SCIPY MINIMIZATION PROCEDURE'''
