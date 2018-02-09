#!/usr/bin/env python
# coding=utf-8

# # # imports # # #
import matplotlib.pyplot as plt
import numpy as np

# input files
filename1 = 'bulk_params.txt'


# import and manage data
with open(filename1,"r") as textfile:
    E_cut_line = next(textfile).split(",")
    E_cut = [float(el) for el in E_cut_line]
    k_points_line = next(textfile).split(",")
    k_points = [float(el) for el in k_points_line]
    energies = np.array([len(E_cut), len(k_points)])
    for line in textfile:
        print line



sigma_Al = np.array(np.loadtxt(filename1))
