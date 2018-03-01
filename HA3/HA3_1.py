#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
from numpy.linalg import inv


def main():

    file_loc = 'HA3/HA3_1.txt'
    # Initial parameters
    ni = 4
    alpha_vec = np.array([0.298073, 1.242567, 5.782948, 38.474970])
    C = np.array([1, 1, 1, 1])

    # Unit conversion
    electron_volt_to_hartree = 1 / 27.211385

    # Calculate constant matrix element
    h = np.zeros((ni, ni))
    S = np.zeros((ni, ni))
    Q = np.zeros((ni, ni, ni, ni))
    for p in range(ni):
        for q in range(ni):
            h[p, q] = calc_h(alpha_vec[p], alpha_vec[q])
            S[p, q] = calc_S(alpha_vec[p], alpha_vec[q])
            for r in range(ni):
                for s in range(ni):
                    Q[p, r, q, s] = calc_Q(alpha_vec[p], alpha_vec[r],
                                           alpha_vec[q], alpha_vec[s])
    C = C_norm(C, S)

    energy_criterion = 10**(-5.) * electron_volt_to_hartree
    E_G_tmp = 1
    E_G = 0
    while (abs(E_G - E_G_tmp) > energy_criterion):
        E_G_tmp = E_G
        F = calc_F(h, Q, C)
        C = calc_C(C, F, S)
        C = C_norm(C, S)
        E_G = calc_ground_state_energy(h, Q, C)
        print(E_G)

    h = 0.005
    rmin = 0
    rmax = 7
    r = np.arange(rmin, rmax, h)

    phi = np.zeros(np.size(r))
    for i in range(len(C)):
        phi = phi - C[i] * np.exp(-alpha_vec[i] * r**2)  # no commnent
    n = phi**2 * 4 * np.pi * r**2
    n = n / np.trapz(n, r)

    with open(file_loc, 'w') as textfile:
        text = str(E_G) + '\n'
        for i in range(len(n)):
            text = text + str(n[i]) + ',' + str(r[i]) + '\n'
        text = text[:-1]
        textfile.write(text)

    with open('HA3/HA3_1_C.txt', 'w') as textfile:
        text = ''
        for val in C:
            text = text + str(val) + '\n'
        textfile.write(text)


def calc_C(C, F, S):
    (E_vec, C_vec) = eig(F, S)
    n_min = np.argmin(E_vec)
    return np.array(C_vec[:, n_min])


def calc_ground_state_energy(h, Q, C, ni=4):
    E_G = 0
    for p in range(ni):
        for q in range(ni):
            E_G += 2 * C[p] * C[q] * h[p, q]
            for r in range(ni):
                for s in range(ni):
                    E_G += Q[p, r, q, s] * C[p] * C[q] * C[r] * C[s]
    return E_G


def calc_F(h, Q, C, ni=4):
    F = np.zeros(np.shape(h))
    for p in range(ni):
        for q in range(ni):
            F[p, q] += h[p, q]
            for r in range(ni):
                for s in range(ni):
                    F[p, q] += Q[p, r, q, s] * C[r] * C[s]
    return F


def calc_Q(alpha_p, alpha_r, alpha_q, alpha_s):
    return 2 * np.pi**(5. / 2.) / ((alpha_p + alpha_q) *
                                   (alpha_r + alpha_s) *
                                   np.sqrt(alpha_p + alpha_q + alpha_r + alpha_s))

# def basis_function(r, alpha):
#     return np.e**(-alpha*r**2)


def calc_S(alpha_p, alpha_q):
    return (np.pi / (alpha_p + alpha_q))**(3. / 2.)


def calc_h(alpha_p, alpha_q):
    return 3 * alpha_q * alpha_p * (np.pi)**(3. / 2.) \
        / ((alpha_p + alpha_q)**(5. / 2.)) - 4 * np.pi / (alpha_p + alpha_q)


def C_norm(C, S, ni=4):
    summ = 0
    for p in range(ni):
        for q in range(ni):
            summ += C[p] * S[p, q] * C[q]
    return C / np.sqrt(summ)


def save_data(name, C, E_G):
    with open('HA3/' + name + '.txt', 'w') as textfile:
        text = 'Ground state energy: ' + str(E_G) + '\n'
        text = text + 'Scaling parameters: '
        for val in C:
            text = text + str(val) + ','
        text = text[:-1]
        textfile.write(text)


if __name__ == '__main__':
    main()
