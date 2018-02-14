import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
def main():
    ''' Finite difference geometry (1D) '''
    r_max = 8 # Maximum radius of position grid in Hartree units
    r_min = 0.0 # Minimum radius of position grid in Hartree units
    n_r = 1000 # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r+2) # Position grid in Hartree units
    r = r[1:-1]
    print(len(r))

    ''' Diffrentiation '''
    delta_r = r[1]-r[0]
    A_dd = create_second_derivative_matrix_1D(r, delta_r)

    ''' Physical constants '''
    Z_helium = 2 # Charge of helium nucleus in hartree units
    Z_hydrogen = 1

    ''' Boundary conditions '''
    U_i = 0
    U_f = 1

    ''' Analytical solutions '''
    V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
    phi_s_H = (1/np.sqrt(np.pi))*np.exp(-r) # Wave function for hydrogen
    f = np.sqrt(4*np.pi)*phi_s_H*r
    n_s_H = phi_s_H*phi_s_H # Electron density hydrogen ground state

    ''' Initial conditions for U '''
    R_U = np.zeros(np.size(r))
    R_U[0] = U_i
    R_U[-1] = U_f

    ''' Initial condition for f '''
    R_f = np.zeros(np.size(r))
    R_f[0] = 0
    R_f[-1] = 0

    ''' Simulation '''
    U = npl.solve(A_dd, -4*np.pi*r*n_s_H - R_U/(delta_r**2))
    A = ((-1/2)*A_dd - np.identity(n_r)*Z_hydrogen/r + R_f*np.identity(n_r)/(2*delta_r**2))


    (energy, wave_functions) = sparse.linalg.eigs(A, which='SM')
    n_min = np.argmin(energy)
    f_min = wave_functions[:, n_min]
    f_min = f_min/np.sqrt(np.trapz(f_min**2, r))

    phi_min = f_min/(np.sqrt(4*np.pi)*r)
    print("dsadsa", trapz(f_min**2, r))
    print("fasfaf", trapz(4*np.pi*r**2*phi_min**2, r))
    phi_min[0] = phi_min[1]/(1-Z_hydrogen*delta_r)
    energy_min = energy[n_min]

    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    ax_potential.plot(r, np.sqrt(phi_min), label="Energy: " + str(2*energy_min))
    ax_potential.plot(r, 1/np.sqrt(np.pi)*np.exp(-r), label="Energy: " + str(2*energy_min))
    ax_potential.set_xlabel('Radius [a.u.]')
    ax_potential.set_ylabel('Probability [a.u]')
    plt.legend()

    with open("simulation_outputs.txt", "w") as textfile:
        textfile.write("Minimum energy: " + str(2*energy_min) + "\n")
        textfile.write("Normalization: " + str(trapz(phi_min**2, r)) + "\n")

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
    print
    return A/(delta_r**2)

if __name__ == '__main__':
    main()
