#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 21:19:47 2019

@author: hao
"""

import numpy as np
import fun_generators as fg
import GLSW

h = np.array([0.0, 1.0, 3.0, 4.0, 4.74])
lenh = len(h)
xs = range(lenh)
q = np.array([0.0, 0.0, 0.0])
ham = GLSW.sw_hamiltonian(q)

for flag in xs:
    field = h[flag]
    fname = 'gs_info/h=' + str(field) + 'T/opt_angles.txt'
    angles = np.load(fname)
    print('h= ', field, '\n')

    for flag1 in range(4):
        theta = angles[flag1*4+2]
        phi   = angles[flag1*4+3]
        Sz = fg.fun_sz(theta, phi)
        Szsq = fg.fun_sz_sq(theta)
        print('sublattice: ',  flag1, '\n')
        print('Sz = ', Sz, 'Sz**2= ', Szsq, '\n')

