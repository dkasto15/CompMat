import numpy as np
import matplotlib.pyplot as plt


def main():
    # Boundary conditions for poisson equation
    U_0 = 0
    U_inf = 0

    # Discretization of 1D grid
    r_max = 10
    n_points = 1000  # Number of points in grid
    r = np.linspace(0, 1, n_points) * r_max  # Coordinate vector
    h = r[1] - r[0]  # Stepsize

    # Differentiation matrix
    T = create_differetiation_matrix(n_points)
    print(T)


def create_differetiation_matrix(n_points):
    # http://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec32.pdf
    T = np.zeros((n_points, n_points))
    for i in range(n_points):
        if i != 0:
            T[i, i - 1] = 1
        if i != n_points - 1:
            T[i, i + 1] = 1
        T[i, i] = -2
    return T


if __name__ == '__main__':
    main()
