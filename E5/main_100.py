import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
def main():
    ''' Finite difference geometry (1D) '''
    r_max = 10 # Maximum radius of position grid in Hartree units
    r_min = 0.0 # Minimum radius of position grid in Hartree units
    n_r = 100 # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r+2) # Position grid in Hartree units
    r = r[1:-1]

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
    f = np.sqrt(4*np.pi)*phi_s_H*r # Assisting wave function
    n_s_H = phi_s_H*phi_s_H # Electron density hydrogen ground state

    R_U = np.zeros(np.size(r))
    R_U[0] = U_i
    R_U[-1] = U_f



    ## LOOP ##
    energy_tmp = 0
    energy_min = 100
    #f = np.sqrt(Z_helium**3*np.exp(-2*Z_helium*r)/np.pi) # Initial guess for wave function

    phi = f/(r*np.sqrt(4*np.pi))
    n_s = 2*(phi)*np.conjugate(phi)


    #V_H = npl.solve(A_dd, -4*np.pi*r*n_s - R_U/(delta_r**2))/r

    R_f = np.zeros(np.size(r))
    R_f[0] = 0
    R_f[-1] = 0

    while(abs(energy_tmp - energy_min) > 0.0001):
        energy_tmp = energy_min

        V_H = npl.solve(A_dd, -4*np.pi*r*n_s - R_U/(delta_r**2))/r

        H = ((-1/2)*A_dd - np.identity(n_r)*Z_helium/r + V_H)

        (energy, wave_functions) = sparse.linalg.eigs(H, 6, which='SM')

        n_min = np.argmin(energy)

        f = wave_functions[:, n_min]

        f = phi/np.sqrt(np.trapz(f**2, r))

        phi = f/(r*np.sqrt(4*np.pi))
        phi[0] = phi[1]/(1-Z_hydrogen*delta_r)
        phi = phi/np.sqrt(trapz(phi*np.conjugate(phi),r))
        n_s = 2*(phi)*np.conjugate(phi)

        energy_min = 2*energy[n_min] - trapz(4*np.pi*r**2*V_H*n_s, r);
        print(energy[n_min])
        #print(energy_min)

    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    print( phi*np.conjugate(phi))
    ax_potential.plot(r, 2*phi*np.conjugate(phi), '*')
    ax_potential.plot(r, 4*r**2*Z_helium**3*np.exp(-2*r*Z_helium))
    ax_potential.plot(r, 4*r**2*(27/16)**3*np.exp(-2*r*Z_helium))
    ax_potential.set_xlabel('Radius [a.u.]')

    with open("simulation_outputs.txt", "w") as textfile:
        textfile.write("Minimum energy: " + str(energy_min) + "\n")
        textfile.write("Normalization: " + str(trapz(phi**2, r)) + "\n")

    plt.show()


def matmul(A, b):
    n = len(b)
    c = np.zeros(np.shape(b))
    for i in range(n):
        for j in range(n):
            c[i] += A[i,j]
    return c

def create_second_derivative_matrix_1D(r, delta_r):
    n = np.size(r)
    A = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i, j] = -2/(delta_r**2)
                if i != 0:
                    A[i, j-1] = 1/(delta_r**2)
                if i != n-1:
                    A[i, j+1] = 1/(delta_r**2)
    return A

if __name__ == '__main__':
    main()
