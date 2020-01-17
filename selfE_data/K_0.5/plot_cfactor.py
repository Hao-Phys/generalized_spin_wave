#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_cfactor.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.16.2020
# Last Modified Date: 01.16.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt('cfactor.txt')
factor = data[:, 0]
colorvec = np.array(['ro', 'g+', 'bv', 'k*', 'c>', 'ms', 'yp', 'k.'])

fig, ax = plt.subplots()

for flag in range(len(factor)):
    for flag1 in range(8):
        plt.plot(factor, abs(data[:, 5+2*flag1]), colorvec[flag1])

plt.show()


