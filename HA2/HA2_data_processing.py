
from ase.io import read
from ase.units import kJ
from ase.utils.eos import EquationOfState
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
<<<<<<< HEAD
import matplotlib.pyplot as plt
=======

# conﬁgs = read('out.txt@:' + str(N_lattice_spacings))  # read 7 conﬁgurations
>>>>>>> 4455cd8a5ce41bb446d8c426a901089040156ea6
conﬁgs = read('HA2/out.txt@0:7')  # read 7 conﬁgurations

# Extract volumes and energies:
volumes = [atoms.get_volume() for atoms in conﬁgs]
energies = [atoms.get_potential_energy() for atoms in conﬁgs]
eos = EquationOfState(volumes, energies)
print(volumes)
print(energies)
plt.plot(volumes)
plt.show()
v0, E_bulk, B = eos.ﬁt()
eos.plot('Al_eos.png')

a_calc = (4 * v0)**(1 / 3.0)  # Is this correct?

N_x = 2
N_y = 2
N_z = 3

surface111 = fcc111('Al', a=a_calc, size=(N_x, N_y, N_z))
surface100 = fcc100('Al', a=a_calc, size=(N_x, N_y, N_z))

surface111.set_calculator(asdasdasd)
surface100.set_calculator(asdasdasd)

cell111 = surface111.get_cell()
area111 = np.linalg.norm(np.cross(cell[0], cell[1]))
sigma111 = (1 / (2 * area111)) * (surface111.get_potential_energy() - N_x * N_y * E_bulk)

cell100 = surface111.get_cell()
area100 = np.linalg.norm(np.cross(cell[0], cell[1]))
sigma100 = (1 / (2 * area100)) * (surface100.get_potential_energy() - N_x * N_y * E_bulk)


atoms = wulff_construction('Al',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1)],
                           energies=[sigma100, sigma111],
                           size=10000,
                           structure='fcc',
                           rounding='below')  # What does this one do?
atoms.center(vacuum=10)
view(atoms)
