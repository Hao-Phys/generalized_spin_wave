#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : main_glsw_disp.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.14.2020
# Last Modified Date: 01.14.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import GLSW
from matplotlib import pyplot as plt


K = - 0.5
inputpath1 = 'cuts/K=' + str(K) + '/path_domain1.dat'
field = cf.field
#inputpath2 = 'cuts/K=' + str(K) + '/path_domain2.dat'
#inputpath3 = 'cuts/K=' + str(K) + '/path_domain3.dat'

domain1 = np.loadtxt(inputpath1)
# domain2 = np.loadtxt(inputpath2)
# domain3 = np.loadtxt(inputpath3)

lenk1 = len(domain1[:, 0])

k = np.arange(lenk1)
omega = np.zeros((8, lenk1))

for flag in range(lenk1):
    
    q1 = domain1[flag, 0]
    q2 = domain1[flag, 1]
    q3 = domain1[flag, 2]
    q = np.array([q1, q2, q3])
    ek, tmp = GLSW.eigensystem(q)
    omega[:, flag] = ek[:8]

fname = 'disp.txt'
np.savetxt(fname, omega)

for band in range(8):
    plt.plot(k, omega[band, :])
plt.show()
