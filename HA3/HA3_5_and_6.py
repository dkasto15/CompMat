import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
import numpy.linalg as npl
import scipy.linalg as spl
from HA3_2 import create_second_derivative_matrix_1D
from HA3_2 import compute_VsH_and_U
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve


def main():
    ''' Run task 4, 5 or 6? '''
    task = 6

    ''' Physical constants '''
    Z_helium = 2  # Charge of helium nucleus in hartree units
    Z_hydrogen = 1  # Charge of hydrogen nucleus in hartree units

    ''' Computation '''
    nbr_of_conv_loops = 1
    Z = Z_helium
    E_vec = np.zeros(nbr_of_conv_loops)
    # eps_vec = np.zeros(nbr_of_conv_loops)
    r_max_vec = np.zeros(nbr_of_conv_loops)
    # phi_vec = np.zeros(nbr_of_conv_loops)

    for j in range(nbr_of_conv_loops):
        ''' Finite difference geometry (1D) '''
        r_max = 7 + j  # Maximum radius of position grid in Hartree units
        r_min = 0  # Minimum radius of position grid in Hartree units
        r_max_vec[j] = r_max
        # n_r = 1000 # Number of elements in position grid
        h = 0.005  # step size
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
        n_s_H = 2*phi_s_H_theor**2  # Initial guess
        n_s_H = n_s_H / trapz(n_s_H * 4 * np.pi * r**2, r) # Normalize

        E_0 = 1
        E_0_old = 0
        counter = 0
        while(abs(E_0 - E_0_old) > 10**(-5) or counter < 3):  # Run at least
              # three iterations to reduce susceptibility to initial values
            counter += 1
            E_0_old = E_0
            # print(trapz(n_s_H * 4 * np.pi * r**2, r)) # Check if normalized
            V_s_H = compute_VsH_and_U(A_dd, r, n_s_H) # Hartree potential
            if task == 4:
                V_H = V_s_H
                V_xc = 0
                eps_xc = 0
            if task == 5:
                V_H = 2 * V_s_H
                V_x, eps_x = calc_xc(n_s_H)[0:2]
                V_xc = V_x
                eps_xc = eps_x
            if task == 6:
                V_H = 2 * V_s_H
                V_x, eps_x, V_c, eps_c, V_xc, eps_xc = calc_xc(n_s_H)
            V = V_H + V_xc
            A = create_Hamiltonian(r, h, Z, V)
            n_s_H, eps, u = compute_n(A, r)
            E_0 = 2 * eps - 2 * trapz(abs(u)**(2.0) * \
                           (0.5 * V_H + V_xc - eps_xc), r) # Ground state energy
        E_vec[j] = E_0 # Energy
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
    ax_pot2.plot(r, 4 * np.pi * r**2 * n_s_H, label='Eigenvalues')
    ax_pot2.set_title('Radial probability distribution for helium electrons when the \n' +
                            'maximum radius in simulation was ' + str(r_max) + ' a.u.')
    ax_pot2.set_xlabel('Radius [a.u.]')
    ax_pot2.set_ylabel('Electron density [a.u.]')

    plt.show()

    ''' Write data to file '''
    print(E_vec)
    with open("eigAndEn.txt", "w") as textfile:
        for j in range(nbr_of_conv_loops):
            textfile.write(str(E_vec[j]) + '\n')

''' Functions '''

def create_Hamiltonian(r, h, Z, V):
    # http://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec32.pdf
    n = np.size(r)
    diag_main = np.ones(n) / (h**2) - Z / r + V
    diag_side = -1 / 2.0 * np.ones(n-1) / (h**2)
    diagonals = [diag_main, diag_side, diag_side]
    A = diags(diagonals, [0, -1, 1])
    return A

def compute_n(A, r):
    (eig, wave) = sparse.linalg.eigs(A, which='SM')  # SM = smallest
    # # magnitude of the eigenvectors
    # (eig, wave) = spl.eig(A) # Find eigenvalues and eigenvectors to Kohn-Sham eq.
    eps_min_ind = np.argmin(eig) # Pick the eigenvalue corresponding to the lowest eigenvalue
    eps = eig[eps_min_ind]
    u = wave[:, eps_min_ind]
    n_s_H = (abs(u) / r)**2 / (4 * np.pi)
    norm = trapz(n_s_H * 4 * np.pi * r**2, r) # Normalization factor
    u = u / np.sqrt(norm)
    n_s_H = n_s_H / norm
    return n_s_H, eps, u

def calc_xc(n):
    ''' Exchange '''
    ex = -(3 / 4.0) * (3 * n / np.pi)**(1 / 3.0)
    dex_dn = -(1 / 4.0) * (3 / np.pi)**(1 / 3.0) * n**(-2 / 3.0)
    V_x = ex + n * dex_dn
    ''' Correlation '''
    A, B, C, D = 0.0311, -0.048, 0.0020, -0.0116
    gamma, beta1, beta2 = -0.1423, 1.0529, 0.3334
    ec = np.zeros(len(n))
    dec_dn = np.zeros(len(n))
    r = (3 / (4 * np.pi * n))**(1 / 3.0)
    drdn = -(1 / 3) * (3 / (np.pi * 4))**(1 / 3) * n**(-4 / 3.0)
    for ind, r_s in enumerate(r):
        if r_s < 1:
            ec[ind] = A * np.log(r_s) + B + C * r_s * np.log(r_s) + D * r_s
            dec_dn[ind] = drdn[ind] * (A / r_s + C * np.log(r_s) + C + D)
        if r_s >= 1:
            ec[ind] = gamma / (1 + beta1 * np.sqrt(r_s) + beta2 * r_s)
            dec_dn[ind] = drdn[ind] * -gamma * (beta1 / (2 * np.sqrt(r_s)) + beta2) \
                        / ((beta1 * np.sqrt(r_s) + beta2 * r_s + 1)**2)
    V_c = ec + n * dec_dn
    ''' Exchange-Correlation '''
    V_xc = V_x + V_c
    e_xc = ex + ec
    return V_x, ex, V_c, ec, V_xc, e_xc

if __name__ == '__main__':
    main()
