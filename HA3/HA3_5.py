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
    nbr_of_conv_loops = 7
    Z = Z_helium
    E_vec = np.zeros(nbr_of_conv_loops)
    eps_vec = np.zeros(nbr_of_conv_loops)
    r_max_vec = np.zeros(nbr_of_conv_loops)
    # phi_vec = np.zeros(nbr_of_conv_loops)

    for j in range(nbr_of_conv_loops):
        ''' Finite difference geometry (1D) '''
        r_max = 3 + 3 * j  # Maximum radius of position grid in Hartree units
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
        V_xc = calc_Vxc_and_epsxc(r)[0]
        eps_xc = calc_Vxc_and_epsxc(r)[1]
        # V_x = 0
        # V_c = 0
        while(abs(E_0 - E_0_old) > 10**(-5) or counter < 3):  # Run at least
              # three iterations to reduce susceptibility to initial values
            counter += 1
            E_0_old = E_0
            phi = u / (np.sqrt(4.0 * np.pi) * r)
            # n_s_H = np.absolute(phi * phi)  # Radial density helium ground state
            # print(trapz(n_s_H * 4 * np.pi * r**2, r)) # Check if normalized
            V_s_H = compute_VsH_and_U(A_dd, r, phi)[0]
            V_H = 2*V_s_H
            A = (-1 / 2.0) * A_dd - I * (Z / r) + V_H + V_xc
            eps, u, E_0 = compute_eps_and_phi(A, r, V_xc, eps_xc, A_dd)  # min eigenvalue & eigenvector
            # print(E_0)
        E_vec[j] = E_0
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
    ax_potential2.set_title('Electron density for helium electons when the \n' +
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

def compute_eps_and_phi(A, r, V_xc, eps_xc, A_dd):
    # (eig, wave) = sparse.linalg.eigs(A, which='SM')  # SM = smallest
    # # magnitude of the eigenvectors
    (eig, wave) = npl.eig(A)
    eig_vec = []
    wave_vec = []
    for ind, e in enumerate(eig): # Remove non-physical solutions (imag. eigenvalues)
        if np.isreal(e):
            eig_vec.append(e)
            wave_vec.append(wave[:, ind])
    eig = np.transpose(eig_vec)
    wave = np.transpose(wave_vec)
    E_0_vec = np.zeros(len(eig))
    for ind, e in enumerate(eig):
        # wave[:,ind] = wave[:,ind] / np.sqrt(np.trapz(np.absolute(wave[:,ind])**2, r))
        phi = wave[:, ind] / (np.sqrt(4.0 * np.pi) * r)
        n_s_H = np.absolute(phi * phi)  # Radial density helium ground state
        norm = np.sqrt(trapz(n_s_H * 4 * np.pi * r**2, r)) # Check if normalized
        wave[:, ind] = wave[:, ind] / norm
        phi = phi / norm
        V_H = compute_VsH_and_U(A_dd, r, phi)[0]
        E_0_vec[ind] = 2 * e - 2 * trapz(np.absolute(wave[:,ind])**(2.0) * (0.5 * V_H + V_xc - eps_xc), r)
    E_min_ind = np.argmin(E_0_vec) # index of lowest energy
    # print(E_min_ind)
    # print(np.argmin(eig))
    eps_min = eig[E_min_ind]
    E = E_0_vec[E_min_ind]
    u = wave[:, E_min_ind]
    # norm = np.sqrt(np.trapz(np.absolute(wave[:, E_min_ind])**2, r))
    # u = wave[:, E_min_ind] / norm # Find the wave function set
    # corresponding to the lowest energy
    # u = u / np.sqrt(np.trapz(u**2, r))  # normalization
    return eps_min, u, E

def calc_Vxc_and_epsxc(r):
    A, B, C, D = 0.0311, -0.048, 0.0020, -0.0116
    gamma, beta1, beta2 = -0.1423, 1.0529, 0.3334
    # ec = []
    # dexc_dn = []
    n = (3 / 4 * np.pi * r**3)
    drdn = -(1 / 3) * (3 / np.pi * 4)**(1 / 3) * n ** (-4 / 3)
    ex = (-(3 / 4) * (3 * n / np.pi)**(1 / 3))
    ec_over_one = gamma / (1 + beta1 * np.sqrt(r) + beta2 * r)
    ec_under_one = A * np.log(r) + B + C * r * np.log(r) + D * r
    dexdn = -(1 / 4) * (3 / np.pi)**(1 / 3) * n**(-2 / 3)
    dec_over_dn = drdn * -gamma * (beta1 / (2 * np.sqrt(r)) + beta2) \
                / ((beta1 * np.sqrt(r) + beta2 * r + 1)**2)
    dec_under_dn = drdn * (A / r + C * np.log(r) + C + D)
    r_val = 0
    counter = 0
    while r_val < 1:
        counter += 1
        r_val = r[counter]
    ec = np.array([])
    ec = np.append(ec, ec_under_one[:counter])
    ec = np.append(ec, ec_over_one[counter:])
    exc = ex + ec
    dexc_dn = np.array([])
    dexc_dn = np.append(dexc_dn, dexdn[:counter] + dec_under_dn[:counter])
    dexc_dn = np.append(dexc_dn, dec_over_dn[counter:] + dexdn[counter:])
    V_xc = exc + n * dexc_dn
    return V_xc, exc

if __name__ == '__main__':
    main()
