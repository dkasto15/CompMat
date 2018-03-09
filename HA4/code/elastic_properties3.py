import matplotlib.pyplot as plt


def main():
    ''' Read in parameters for EAM calc '''
    with open('HA4/results/fit_potential_output.txt', 'r') as textfile:
        line = next(textfile)
        line = line.split(',')

        A = float(line[0])
        lmbd = float(line[1])
        D = float(line[2])
        mu2 = float(line[3])

    path_BCC = 'HA4/downloads/Al-EV-curves/BCC.dat'
    V_BCC, E_BCC, scale_factor = read_data(path_BCC)

    path_DI = 'HA4/downloads/Al-EV-curves/DIAMOND.dat'
    V_DI, E_DI, scale_factor = read_data(path_DI)

    path_FCC = 'HA4/downloads/Al-EV-curves/FCC.dat'
    V_FCC, E_FCC, scale_factor = read_data(path_FCC)

    path_SC = 'HA4/downloads/Al-EV-curves/SC.dat'
    V_SC, E_SC, scale_factor = read_data(path_SC)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(V_BCC, E_BCC)
    ax.plot(V_DI, E_DI)
    ax.plot(V_FCC, E_FCC)
    ax.plot(V_SC, E_SC)
    plt.show()


def read_data(path):
    E = []
    V = []
    scale_factor = []
    with open(path, 'r') as textfile:
        line = next(textfile)
        a0_BCC = float(line.split(':')[1])
        next(textfile)
        next(textfile)
        next(textfile)
        next(textfile)
        for line in textfile:
            line = line.split()
            V.append(float(line[0]))
            E.append(float(line[1]))
            scale_factor.append(float(line[2]))

    return V, E, scale_factor


if __name__ == '__main__':
    main()
