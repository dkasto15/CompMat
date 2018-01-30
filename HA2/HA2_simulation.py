
#!/usr/bin/env python
# # imports # #
import numpy as np
from gpaw import GPAW, Mixer, PW
from ase.build import bulk
from ase.io import read
from ase.units import kJ
from ase.utils.eos import EquationOfState
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view

a_al = 4.05
al = bulk('Al', 'fcc', a=a_al, cubic=False)

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)
N_lattice_spacings = 7
E_al = 84.67567  # ionization energy for hardest bound core electron
E_cut = [50, E_al, 100, 200, 300, 400, 500, 600, 700, 800]  # cut-off energy

for energy in [500]:
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=mixer,
                kpts=(12, 12, 12),  # k-point grid
                txt='out.txt')  # name of GPAW output text file

al.set_calculator(calc)
cell_0 = al.cell
for eps in np.linspace(-0.02,	0.02,	N_lattice_spacings):
    al.cell = (1 + eps) * cell_0
    al.get_potential_energy()


conﬁgs = read('out.txt@:' + str(N_lattice_spacings))  # read 7 conﬁgurations

# Extract volumes and energies:
volumes = [atoms.get_volume() for atoms in conﬁgs]
energies = [atoms.get_potential_energy() for atoms in conﬁgs]
eos = EquationOfState(volumes, energies)
v0, E_bulk, B = eos.ﬁt()
eos.plot('Al_eos.png')

a_calc = (4 * v0)**(1 / 3.0)  # Is this correct?

N_x = 2
N_y = 2
N_z = 3

surface111 = fcc111('Al', a=a_calc, size=(N_x, N_y, N_z))
surface100 = fcc100('Al', a=a_calc, size=(N_x, N_y, N_z))

surface111.set_calculator(calc)
surface100.set_calculator(calc)

cell111 = surface111.get_cell()
area111 = np.linalg.norm(np.cross(cell[0], cell[1]))
sigma111 = (1 / (2 * area111)) * (surface111.get_potential_energy() - N_x * N_y * E_bulk)

cell100 = surface111.get_cell()
area100 = np.linalg.norm(np.cross(cell[0], cell[1]))
sigma100 = (1 / (2 * area100)) * (surface100.get_potential_energy() - N_x * N_y * E_bulk)


al_construction = wulff_construction('Al',
                                     surfaces=[(1, 0, 0),
                                               (1, 1, 1)],
                                     energies=[sigma100, sigma111],
                                     size=10000,
                                     structure='fcc',
                                     rounding='below')  # Vad gör denna?
al_construction.center(vacuum=10)


view(atoms)
