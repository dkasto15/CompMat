import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags, csc_matrix


def main():
    ''' Finite difference geometry (1D) '''
    r_max = 10  # Maximum radius of position grid in Hartree units
    r_min = 0.01  # Minimum radius of position grid in Hartree units
    n_r = 1000  # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r)  # Position grid in Hartree units

    ''' Differentiation '''
    h = r[1] - r[0]  # Step size
    # print(h)
    A = create_second_derivative_matrix_1D(r, h)

    ''' Physical constants '''
    Z_helium = 2  # Charge of helium nucleus in hartree units
    Z_hydrogen = 1  # Charge of hydrogen nucleus in hartree units

    ''' Boundary conditions for poisson equation '''
    # U_0 = 0
    # U_inf = 0

    ''' Analytical solutions '''
    V_hartree = 1 / r - (1 + 1 / r) * np.exp(-2 * r) # The hartree potential for hydrogen
    n_s_H = ((1 / np.sqrt(np.pi)) * np.exp(-r))**2  # Electron density hydrogen ground state

    ''' Computation '''
    V_s_H = compute_VsH_and_U(A, r, n_s_H)

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    ax_potential.plot(r, V_s_H, label='Calculated V$_{H}$')
    ax_potential.plot(r, V_hartree, '--', label='Theoretical V$_{H}$')
    ax_potential.set_xlabel('Radius [a.u.]')
    ax_potential.set_ylabel('Hartree potential [a.u.]')
    ax_potential.set_title('Difference between simulated and analytical potential')
    ax_potential.legend(loc=1)

    plt.show()

''' Functions '''


def compute_VsH_and_U(A, r, n_s_H):
    # print(trapz(n_s_H*4*np.pi*r**2, r)) # Check if normalized
    U_0 = spsolve(A, -4 * np.pi * r * n_s_H)
    U = U_0 + r / r[-1]
    V_s_H = U / r
    return V_s_H


def create_second_derivative_matrix_1D(r, h):
    # http://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec32.pdf
    n = np.size(r)
    diag_main = np.ones(n) * -2 / (h**2)
    diag_side = np.ones(n-1) / (h**2)
    diagonals = [diag_main, diag_side, diag_side]
    A = csc_matrix(diags(diagonals, [0, -1, 1]))
    return A


if __name__ == '__main__':
    main()
