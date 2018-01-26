from ase.io import read
from ase.units import kJ
from ase.utils.eos import EquationOfState
from ase.build import fcc111, fcc100
from ase.cluster.wulff import wulff_construction
from ase.visualize import view

conﬁgs = read('out.txt@0:7')  # read 5 conﬁguratons
# Extract volumes and energies:
volumes = [atoms.get_volume() for atoms in conﬁgs]
energies = [atoms.get_potential_energy() for atoms in conﬁgs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.ﬁt()
eos.plot('Al_eos.png')

a_calc = (4 * v0)**(1 / 3.0)

surface111 = fcc111('Al', a=a_calc, size=(2, 2, 3))
surface100 = fcc100('Al', a=a_calc, size=(2, 2, 3))


atoms = wulff_construction('Al',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1)],
                           energies=[0.6, 0.6],
                           size=10000,
                           structure='fcc',
                           rounding='below')
atoms.center(vacuum=10)
view(atoms)
