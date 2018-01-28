
#!/usr/bin/env python

# # imports # #
from gpaw import GPAW, Mixer
from ase.build import bulk

a_al = 4.05
al = bulk('Al', 'fcc', a=a_al, cubic=False)

mixer = Mixer(beta=0.1,	nmaxold=5,	weight=50.0)

E_al = 84.67567
for energy in [E_al, 100, 200, 300, 400, 500, 600, 700, 800]:
    calc = GPAW(mode=PW(energy),  # use the LCAO basis mode
                h=0.18,  # grid spacing
                xc='PBE',  # XC-functional
                mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0)
                kpts=(12, 12, 12),  # k-point grid
                txt='out.txt')  # name of GPAW output text file
al.set_calculator(calc)

cell_0 = al.cell
for eps in np.linspace(-0.02,	0.02,	7):
    al.cell = (1 + eps) * cell0
    al.get_potential_energy()
