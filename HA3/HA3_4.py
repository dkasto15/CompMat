import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
import numpy.linalg as npl
from HA3_2 import create_second_derivative_matrix_1D
from HA3_2 import compute_VsH_and_U


def main():
    ''' Physical constants '''
    Z_helium = 2  # Charge of helium nucleus in hartree units
    Z_hydrogen = 1  # Charge of hydrogen nucleus in hartree units

    ''' Computation '''
    nbr_of_conv_loops = 1
    Z = Z_helium
    E_vec = np.zeros(nbr_of_conv_loops)
    eps_vec = np.zeros(nbr_of_conv_loops)
    r_max_vec = np.zeros(nbr_of_conv_loops)
    # phi_vec = np.zeros(nbr_of_conv_loops)

    for j in range(nbr_of_conv_loops):
        ''' Finite difference geometry (1D) '''
        r_max = 10 + 2 * j  # Maximum radius of position grid in Hartree units
        r_min = 0  # Minimum radius of position grid in Hartree units
        r_max_vec[j] = r_max
        # n_r = 1000 # Number of elements in position grid
        h = 0.1  # step size, based on that it was sufficient for task 2
        r = np.arange(r_min, r_max, h)  # constant step size
        # r = np.linspace(r_min, r_max, n_r+1) # Position grid in Hartree units
        r = r[1:]  # Remove singularity in r=0
        I = np.identity(len(r))

        ''' Differentiation '''
        A_dd = create_second_derivative_matrix_1D(r, h)

        ''' Analytical solutions (and initial guesses) '''
        # V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
        phi_s_H_theor = (1 / np.sqrt(np.pi)) * np.exp(-r)  # Wave function for hydrogen in hartree
        u = phi_s_H_theor  # initial guess

        E_0 = 1
        E_0_old = 0
        counter = 0
        V_xc = 0
        eps_xc = 0
        V_x = 0
        V_c = 0
        while(abs(E_0 - E_0_old) > 10**(-5) or counter < 3):  # Run at least
              # three iterations to reduce susceptibility to initial values
            counter += 1
            E_0_old = E_0
            phi = u/(np.sqrt(4*np.pi)*r)
            V_s_H = compute_VsH_and_U(A_dd, r, phi)[0]
            V_H = V_s_H
            A = (-1 / 2.0) * A_dd - I * (Z / r) + V_H + V_x + V_c
            eps, u = compute_eps_and_phi(A, r)  # min eigenvalue & eigenvector
            E_0 = 2 * eps - 2 * trapz(abs(u)**(2.0) * (0.5 * V_H + V_xc - eps_xc), r)
            # print(E_0)

        E_vec[j] = np.real(E_0)  # "real" to suppress annoying complex warning
        eps_vec[j] = eps
        print('Done: r_max = ' + str(r_max))
        # phi_vec[j] = phi_s_H

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    #ax_potential.plot(r_max_vec, eps_vec, label='Eigenvalues')
    ax_potential.plot(r_max_vec, E_vec, '--', label='Energies')
    ax_potential.set_xlabel('Max. radius [a.u.]')
    ax_potential.set_ylabel('Converged energy value [a.u]')
    ax_potential.set_title('Energy convergence when increasing maximum radius')
    # ax_potential.legend(loc=2)
    plt.savefig('eigAndEn.eps')
    plt.savefig('eigAndEn.png')

    fig_2 = plt.figure()
    ax_potential2 = fig_2.add_subplot(111)
    ax_potential2.plot(r, 2 * 4 * np.pi * r**2 * u**2, label='Eigenvalues')
    ax_potential2.set_title('Electron densty for helium electons when the \n' +
                            'maximum radius in simulation was ' + str(r_max) + ' a.u.')
    ax_potential2.set_xlabel('Radius [a.u.]')
    ax_potential2.set_ylabel('Electron densinsity [a.u.]')

    plt.show()

    ''' Write data to file '''
    print(eps_vec)
    print(E_vec)
    with open("eigAndEn.txt", "w") as textfile:
        for j in range(nbr_of_conv_loops):
            textfile.write(str(eps_vec[j]) + '\t')
            textfile.write(str(E_vec[j]) + '\n')

''' Functions '''

def compute_eps_and_phi(A, r):
    (eig, wave_functions) = sparse.linalg.eigs(A, which='SM')  # SM = smallest
    # magnitude of the eigenvectors
    # (eig, wave_functions) = npl.eig(A)
    eig = [np.real(i) for i in eig]
    eps_min_ind = np.argmin(eig) # index of lowest eigenvalue (epsilon) in the vector
    eps_min = eig[eps_min_ind]
    u = wave_functions[:, eps_min_ind] # Find the wave function set
    # corresponding to the lowest energy
    # u = u / np.sqrt(np.trapz(u**2, r))  # normalization
    return eps_min, u


def calc_Vxc(r):
    n = (3 / 4 * np.pi * r**3)
    drdn = -(1 / 3) * (3 / np.pi * 4)**(1 / 3) * n ** (-4 / 3)

    ex = -(3 / 4) * (3 * n / np.pi)**(1 / 3)

    dexdn = -(3 / 4) * (3 / np.pi)**(1 / 3) * n**(-1 / 3)

    ec_more_one = gamma / (1 + beta1 * np.sqrt(r) + beta1 * r)
    dec_more_onedr = -gamma * (beta1 / (2 * np.sqrt(r)) + beta2) / (beta1 * np.sqrt())
    ec_less_one = A * np.log(r) + B + C * r * np.log(r) + D * r

    # ex = -(3 / 4) * (3 * 3 / 4 * np.pi)*


if __name__ == '__main__':
    main()
