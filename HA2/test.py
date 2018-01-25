#!/usr/bin/env python
from ase.cluster.wulff import wulff_construction
from ase.visualize import view

atoms = wulff_construction('Al',
                          surfaces=[(1, 0, 0),
                                    (1, 1, 1),
                                    (1, 1, 0)],
                           energies=[0.6, 0.6, 0.6],
                           size=10000,
                           structure='fcc',
                           rounding='below')
atoms.center(vacuum=10)
view(atoms) 
