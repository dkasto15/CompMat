import numpy as np


def main():
    u_0 = 0
    u_max = 0

    # Discretization of 1D grid
    r_max = 10
    n_points = 1000  # Number of points in grid
    r = np.linspace(0, 1, n_points) * r_max  # Coordinate vector
    h = r[1] - r[0]  # Stepsize


if __name__ == '__main__':
    main()
