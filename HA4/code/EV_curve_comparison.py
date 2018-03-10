import matplotlib.pyplot as plt
from ase.build import bulk
from eam_calculator import get_calc


def main():
    ''' Read in parameters for EAM calc '''
    with open('HA4/results/fit_potential_output_full.txt', 'r') as textfile:
        line = next(textfile)
        line = line.split(',')

        A = float(line[0])
        lmbd = float(line[1])
        D = float(line[2])
        mu2 = float(line[3])

    path_BCC = 'HA4/downloads/Al-EV-curves/BCC.dat'
    a0_BCC, V_BCC, E_BCC, scale_factors = read_data(path_BCC)

    path_DI = 'HA4/downloads/Al-EV-curves/DIAMOND.dat'
    a0_DI, V_DI, E_DI, scale_factors = read_data(path_DI)

    path_FCC = 'HA4/downloads/Al-EV-curves/FCC.dat'
    a0_FCC, V_FCC, E_FCC, scale_factors = read_data(path_FCC)
    a0_FCC_experimental = 4.032
    E0_FCC_experimental = -3.36
    path_SC = 'HA4/downloads/Al-EV-curves/SC.dat'
    a0_SC, V_SC, E_SC, scale_factors = read_data(path_SC)

    al_bulk_BCC = bulk('Al', 'bcc', a=a0_BCC)
    al_bulk_DI = bulk('Al', 'diamond', a=a0_DI)
    al_bulk_FCC = bulk('Al', 'fcc', a=a0_FCC)
    al_bulk_SC = bulk('Al', 'sc', a=a0_SC)

    calc = get_calc((A, lmbd, D, mu2))

    al_bulk_BCC.set_calculator(calc)
    al_bulk_DI.set_calculator(calc)
    al_bulk_FCC.set_calculator(calc)
    al_bulk_SC.set_calculator(calc)

    V_BCC_sim, V_DI_sim, V_FCC_sim, V_SC_sim = [], [], [], []
    E_BCC_sim, E_DI_sim, E_FCC_sim, E_SC_sim = [], [], [], []

    for scale_factor in scale_factors:
        cell_0 = al_bulk_BCC.get_cell()
        al_bulk_BCC.set_cell(cell_0 * scale_factor)
        E_BCC_sim.append(al_bulk_BCC.get_potential_energy())
        V_BCC_sim.append(2 * al_bulk_BCC.get_volume())
        al_bulk_BCC.set_cell(cell_0)

        cell_0 = al_bulk_DI.get_cell()
        al_bulk_DI.set_cell(cell_0 * scale_factor)
        E_DI_sim.append(al_bulk_DI.get_potential_energy())
        V_DI_sim.append(4 * al_bulk_DI.get_volume())
        al_bulk_DI.set_cell(cell_0)

        cell_0 = al_bulk_FCC.get_cell()
        al_bulk_FCC.set_cell(cell_0 * scale_factor)
        E_FCC_sim.append(al_bulk_FCC.get_potential_energy())
        V_FCC_sim.append(4 * al_bulk_FCC.get_volume())
        al_bulk_FCC.set_cell(cell_0)

        cell_0 = al_bulk_SC.get_cell()
        al_bulk_SC.set_cell(cell_0 * scale_factor)
        E_SC_sim.append(al_bulk_SC.get_potential_energy())
        V_SC_sim.append(1 * al_bulk_SC.get_volume())
        al_bulk_SC.set_cell(cell_0)

    fig = plt.figure(1, figsize=(7.1, 7.1))
    ax_BCC = fig.add_subplot(221)
    ax_BCC.plot(V_BCC, E_BCC, label='DFT', color='blue')
    ax_BCC.plot(V_BCC_sim, E_BCC_sim, label='EAM', color='red')
    ax_BCC.set_title('BCC unit cell')
    ax_BCC.minorticks_on()
    ax_BCC.grid(True, which='minor', linestyle='--')
    ax_BCC.grid(True, which='major', linestyle='-')
    ax_BCC.set_xlim(V_BCC[0], V_BCC[-1])
    ax_BCC.legend()

    ax_DI = fig.add_subplot(222)
    ax_DI.plot(V_DI, E_DI, label='DFT', color='blue')
    ax_DI.plot(V_DI_sim, E_DI_sim, label='EAM', color='red')
    ax_DI.set_title('Diamond unit cell')
    ax_DI.minorticks_on()
    ax_DI.grid(True, which='minor', linestyle='--')
    ax_DI.grid(True, which='major', linestyle='-')
    ax_DI.set_xlim(V_DI[0], V_DI[-1])
    ax_DI.legend()

    ax_FCC = fig.add_subplot(223)
    ax_FCC.plot(V_FCC, E_FCC, label='DFT', color='blue')
    ax_FCC.plot(V_FCC_sim, E_FCC_sim, label='EAM', color='red')
    ax_FCC.plot(a0_FCC_experimental**3, E0_FCC_experimental, marker='*', color='black')
    ax_FCC.set_title('FCC unit cell')
    ax_FCC.minorticks_on()
    ax_FCC.grid(True, which='minor', linestyle='--')
    ax_FCC.grid(True, which='major', linestyle='-')
    ax_FCC.set_xlim(V_FCC[0], V_FCC[-1])
    ax_FCC.legend()

    ax_SC = fig.add_subplot(224)
    ax_SC.plot(V_SC, E_SC, label='DFT', color='blue')
    ax_SC.plot(V_SC_sim, E_SC_sim, label='EAM', color='red')
    ax_SC.set_title('SC unit cell')
    ax_SC.minorticks_on()
    ax_SC.grid(True, which='minor', linestyle='--')
    ax_SC.grid(True, which='major', linestyle='-')
    ax_SC.set_xlim(V_SC[0], V_SC[-1])
    ax_SC.legend()

    fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=None, hspace=0.25)
    fig.text(0.5, 0.03, 'Volume [Å^3]', ha='center', fontsize=16)
    fig.text(0.03, 0.5, 'Energy [eV]', va='center', rotation='vertical', fontsize=16)

    fig.savefig("HA4/results/EV_curves.png", bbox_inches='tight')

    with open('HA4/results/EV_FCC.txt', 'w') as textfile:
        text = 'Volume, Enegy DFT, Energy EAM ' + '\n'
        for i in range(len(V_FCC)):
            text = text + str(V_FCC[i]) + ',' + str(E_FCC[i]) + ',' + str(E_FCC_sim[i]) + '\n'
        text = text[:-1]
        textfile.write(text)

    plt.show()


def read_data(path):
    E = []
    V = []
    scale_factor = []
    with open(path, 'r') as textfile:
        line = next(textfile)
        a0 = float(line.split(':')[1])
        next(textfile)
        next(textfile)
        next(textfile)
        next(textfile)
        for line in textfile:
            line = line.split()
            V.append(float(line[0]))
            E.append(float(line[1]))
            scale_factor.append(float(line[2]))

    return a0, V, E, scale_factor


if __name__ == '__main__':
    main()
