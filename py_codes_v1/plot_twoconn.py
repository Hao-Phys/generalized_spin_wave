#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_twoconn.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.11.2020
# Last Modified Date: 01.12.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
from matplotlib import pyplot as plt
import sys
import GLSW

"""
run this program as: python3 plot_twoconn.py <input..> <min..> <max..>
"""
fname2 = sys.argv[2]
fname3 = sys.argv[3]
num_sub = 4

data_min = np.loadtxt(fname2)
data_max = np.loadtxt(fname3)

lenk = len(data_min[:, 0])
k = np.arange(lenk)
ek = np.zeros((lenk, 2*num_sub))

for flag in range(lenk):

    k1 = data_min[flag, 1]
    k2 = data_min[flag, 2]
    k3 = data_min[flag, 3]

    kk = np.array([k1, k2, k3])
    omegak, tmp = GLSW.eigensystem(kk)
    ek[flag, :] = omegak[:2*num_sub]
    

fig, ax = plt.subplots()

for band in range(2*num_sub):
    plt.plot(k, ek[:, band], 'k-')

for flag1 in range((2*num_sub)**2):
    plt.plot(k, data_min[:, 4+flag1], 'r-.')
    plt.plot(k, data_max[:, 4+flag1], 'b-.')

for band1 in range(2*num_sub):
    for band2 in range(2*num_sub):

        tmin = data_min[:, 4+8*band1+band2]
        tmax = data_max[:, 4+8*band1+band2]

        if (band1 in range(num_sub)) and (band2 in range(num_sub)):
            plt.fill_between(k, tmin, tmax, 
                    facecolor='lavender')

        elif (band1 in range(num_sub, 2*num_sub)) and (band2 in range(num_sub)):
            plt.fill_between(k, tmin, tmax, 
                    facecolor='lightgreen')

        elif (band1 in range(num_sub, 2*num_sub)) and (band2 in range(num_sub)):
            plt.fill_between(k, tmin, tmax, 
                    facecolor='lightblue')

        else:
            plt.fill_between(k, tmin, tmax, 
                    facecolor='wheat')


plt.xticks([0, lenk/4-1, lenk/2-1, 3*lenk/4-1, lenk-1], \
           [-2.0, -1.0, -0.0, 1.0, 2.0]) 
ax.legend()
ax.text(50, 11.0, r'$K=0.5$', fontsize=14)
plt.xlim(0, lenk-1)
plt.ylim(0, 12)
plt.xlabel(r'$H$')
plt.ylabel(r'$\omega$ meV')
plt.show()



