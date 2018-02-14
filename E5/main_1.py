import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
def main():
    ''' Finite difference geometry (1D) '''
    r_max = 5 # Maximum radius of position grid in Hartree units
    r_min = 0.01 # Minimum radius of position grid in Hartree units
    n_r = 1000 # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r) # Position grid in Hartree units

    ''' Diffrentiation '''
    delta_r = r[1]-r[0]
    A = create_second_derivative_matrix_1D(r, delta_r)

    ''' Physical constants '''
    Z_helium = 2 # Charge of helium nucleus in hartree units
    Z_hydrogen = 1

    ''' Boundary conditions '''
    U_i = 0
    U_f = 1

    ''' Analytical solutions '''
    V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
    phi_s_H = (1/np.sqrt(np.pi))*np.exp(-r) # Wave function for hydrogen
    n_s_H = phi_s_H*phi_s_H # Electron density hydrogen ground state
    print(trapz(n_s_H, r))
    R = np.zeros(np.size(r))
    R[0] = U_i
    R[-1] = U_f
    U = npl.solve(A, -4*np.pi*r*n_s_H - R/(delta_r**2))

    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    ax_potential.plot(r, U)
    ax_potential.plot(r, V_hartree*r, '--')
    ax_potential.set_xlabel('Radius [a.u.]')


    plt.show()

def create_second_derivative_matrix_1D(r, delta_r):
    n = np.size(r)
    A = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i, j] = -2
                if i != 0:
                    A[i, j-1] = 1
                if i != n-1:
                    A[i, j+1] = 1
    return A/(delta_r**2)

if __name__ == '__main__':
    main()
