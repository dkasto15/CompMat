
#!/usr/bin/env python

# # imports # #
from gpaw import GPAW
from ase.build import bulk

al = bulk('Al','fcc',a=4.05,cubic=False)

calc = GPAW(mode='lcao',        # use the LCAO basis mode
            basis='dzp',        # use double zeta polarized basis set
            h=0.18,             # grid spacing
            xc='PBE',           # XC-functional
            kpts=(12,12,12),    # k-point grid
            txt='out.txt')      # name of GPAW output text file

al.set_calculator(calc)
al.get_potential_energy()
