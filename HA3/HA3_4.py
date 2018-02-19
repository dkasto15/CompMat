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
    nbr_of_conv_loops = 5
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
        h = 0.2  # step size
        r = np.arange(r_min, r_max, h)  # constant step size
        # r = np.linspace(r_min, r_max, n_r+1) # Position grid in Hartree units
        r = r[1:]  # Remove singularity in r=0
        I = np.identity(len(r))

        ''' Differentiation '''
        A_dd = create_second_derivative_matrix_1D(r, h)

        ''' Analytical solutions (and initial guesses) '''
        # V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
        phi_s_H_theor = (1 / np.sqrt(np.pi)) * np.exp(-r)  # Radial wave function
                                                           # for hydrogen in hartree
        u = phi_s_H_theor  # Initial guess

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
            phi = u / (np.sqrt(4.0 * np.pi) * r) # The wave function
            # n_s_H = np.absolute(phi * phi)  # Radial density helium ground state
            # print(trapz(n_s_H * 4 * np.pi * r**2, r)) # Check if normalized
            V_s_H = compute_VsH_and_U(A_dd, r, phi)[0] # Electro-static potential
            V_H = V_s_H # Hartree potential, approximated as V_s_H
            A = (-1 / 2.0) * A_dd - I * (Z / r) + V_H + V_x + V_c # Kohn-Sham matrix
            eps, u, E_0 = compute_eps_and_phi(A, r, V_xc, eps_xc, A_dd) # min
            # eigenvalue, eigenvector & energy
            # print(E_0)
        E_vec[j] = E_0 # Energy
        eps_vec[j] = eps # Eigenvalue
        print('Done: r_max = ' + str(r_max))
        # phi_vec[j] = phi_s_H

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_pot = fig_1.add_subplot(111)
    #ax_pot.plot(r_max_vec, eps_vec, label='Eigenvalues')
    ax_pot.plot(r_max_vec, E_vec, '--', label='Energies') # Plot E as func. of r_max
    ax_pot.set_xlabel('Max. radius [a.u.]')
    ax_pot.set_ylabel('Converged energy value [a.u]')
    ax_pot.set_title('Energy convergence when increasing maximum radius')
    # ax_pot.legend(loc=2)
    plt.savefig('eigAndEn.eps')
    plt.savefig('eigAndEn.png')

    fig_2 = plt.figure()
    ax_pot2 = fig_2.add_subplot(111)
    ax_pot2.plot(r, 4 * np.pi * r**2 * np.absolute(phi)**2, label='Eigenvalues')
    ax_pot2.set_title('Radial probability distribution for helium electrons when the \n' +
                            'maximum radius in simulation was ' + str(r_max) + ' a.u.')
    ax_pot2.set_xlabel('Radius [a.u.]')
    ax_pot2.set_ylabel('Electron density [a.u.]')

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
    (eig, wave) = npl.eig(A) # Find eigenvalues and eigenvectors to Kohn-Sham eq.
    eig_vec = []
    wave_vec = []
    for ind, e in enumerate(eig): # Remove non-physical solutions (imag. eigenvalues)
        if np.isreal(e):
            eig_vec.append(e)
            wave_vec.append(wave[:, ind])
    eig = np.transpose(eig_vec)
    wave = np.transpose(wave_vec)
    E_0_vec = np.zeros(len(eig))
    for ind, e in enumerate(eig): # Pick the solution with the lowest energy
        phi = wave[:, ind] / (np.sqrt(4.0 * np.pi) * r)
        n_s_H = np.absolute(phi * phi)  # Electron density helium ground state / 2
        norm = np.sqrt(trapz(n_s_H * 4 * np.pi * r**2, r)) # Normalization factor
        wave[:, ind] = wave[:, ind] / norm # Normalize eigenfunction (u)
        phi = phi / norm # Normalize wave function
        V_H = compute_VsH_and_U(A_dd, r, phi)[0] # Hartree potential
        E_0_vec[ind] = 2 * e - 2 * trapz(np.absolute(wave[:,ind])**(2.0) * \
                       (0.5 * V_H + V_xc - eps_xc), r) # Ground state energy
    E_min_ind = np.argmin(E_0_vec) # index of lowest energy
    # print(E_min_ind) # Check if E calc above is necessary or if eig would yield same index
    # print(np.argmin(eig)) # Check if E calc above is necessary or if eig would yield same index
    eps_min = eig[E_min_ind] # Pick the eigenvalue corresponding to the lowest energy
    E = E_0_vec[E_min_ind] # Pick the lowest energy
    u = wave[:, E_min_ind] # Pick the wave function set corresponding to the lowest energy
    return eps_min, u, E

def calc_Vxc_and_epsxc(r): # To be used in task 5
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
