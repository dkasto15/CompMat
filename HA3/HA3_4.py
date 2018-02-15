import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
from HA3_2 import create_second_derivative_matrix_1D
from HA3_2 import compute_VsH_and_U
def main():
    ''' Physical constants '''
    Z_helium = 2 # Charge of helium nucleus in hartree units
    Z_hydrogen = 1 # Charge of hydrogen nucleus in hartree units

    ''' Computation '''
    nbr_of_conv_loops = 8
    Z = Z_hydrogen
    E_vec = np.zeros(nbr_of_conv_loops)
    eps_vec = np.zeros(nbr_of_conv_loops)
    r_max_vec = np.zeros(nbr_of_conv_loops)
    # phi_vec = np.zeros(nbr_of_conv_loops)

    for j in range(nbr_of_conv_loops):
        ''' Finite difference geometry (1D) '''
        r_max = 5+3*j # Maximum radius of position grid in Hartree units
        r_min = 0 # Minimum radius of position grid in Hartree units
        r_max_vec[j] = r_max
        # n_r = 1000 # Number of elements in position grid
        h = 0.1 # step size, based on that it was sufficient for task 2
        r = np.arange(r_min, r_max, h) # constant step size
        # r = np.linspace(r_min, r_max, n_r+1) # Position grid in Hartree units
        r = r[1:] # Remove singularity in r=0
        I = np.identity(len(r))

        ''' Differentiation '''
        A_dd = create_second_derivative_matrix_1D(r, h)

        ''' Analytical solutions (and initial guesses) '''
        # V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
        phi_s_H_theor = (1/np.sqrt(np.pi))*np.exp(-r) # Wave function for hydrogen in hartree
        phi_s_H = phi_s_H_theor # initial guess

        eps = 1
        eps_old = 0
        counter = 0
        while(abs(eps - eps_old) > 10**(-5) or counter < 3): # Run at least
              # three iterations to reduce susceptibility to initial values
            counter += 1
            eps_old = eps
            V_s_H = compute_VsH_and_U(A_dd, r, phi_s_H)[0]
            V_H = V_s_H
            V_x = 0
            V_c = 0
            A = (-1/2.0)*A_dd - I*(Z/r) - V_H - V_x - V_c
            eps, phi_s_H = compute_eps_and_phi(A, r) # min eigenvalue & eigenvector

        V_xc = 0
        eps_xc = 0
        E_0 = 2*eps - 2*trapz(phi_s_H**(2)*(0.5*V_H + V_xc - eps_xc), r)
        E_vec[j] = np.real(E_0) # "real" to suppress annoying complex warning
        eps_vec[j] = eps
        print('Done: r_max = ' + str(r_max))
        # phi_vec[j] = phi_s_H

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    ax_potential.plot(r_max_vec, eps_vec, label='Eigenvalues')
    # ax_potential.plot(r_max_vec, E_vec, '--', label='Energies')
    ax_potential.set_xlabel('Max. radius [atomic units]')
    ax_potential.set_ylabel('Converged value [a.u]')
    ax_potential.legend(loc=1)
    plt.savefig('eigAndEn.eps')
    plt.savefig('eigAndEn.png')
    # plt.show()

    ''' Write data to file '''
    print(eps_vec)
    print(E_vec)
    with open("eigAndEn.txt", "w") as textfile:
        for j in range(nbr_of_conv_loops):
            textfile.write(str(eps_vec[j]) + '\t')
            textfile.write(str(E_vec[j]) + '\n')

''' Functions '''
def compute_eps_and_phi(A, r):
    (eig, wave_functions) = sparse.linalg.eigs(A, which='SM') # SM = smallest
    # magnitude of the eigenvectors
    eig = [np.real(i) for i in eig]
    eps_min_ind = np.argmin(eig) # index of lowest energy in the energy vector
    eps_min = eig[eps_min_ind]
    phi_min = wave_functions[:, eps_min_ind] # Find the wave function set
    # corresponding to the lowest energy
    phi_s_H = phi_min/np.sqrt(np.trapz(phi_min**2, r)) # normalization
    return eps_min, phi_s_H

if __name__ == '__main__':
    main()
