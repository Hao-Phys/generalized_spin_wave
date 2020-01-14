#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 15:34:11 2019

@author: Hao

""
To run this program:
    python3 main_glsw_sqw.py <input.txt>,
    where <input.txt> is the input file for cofig
"""

import numpy as np
import cofig as cf
import GLSW
# from matplotlib import pyplot as plt
# from matplotlib import cm

K = 0.5
inputpath1 = 'cuts/K=' + str(K) + '/path_domain1.dat'
field = cf.field
#inputpath2 = 'cuts/K=' + str(K) + '/path_domain2.dat'
#inputpath3 = 'cuts/K=' + str(K) + '/path_domain3.dat'

domain1 = np.loadtxt(inputpath1)
# domain2 = np.loadtxt(inputpath2)
# domain3 = np.loadtxt(inputpath3)

lenk1 = len(domain1[:, 0])

# the energy range
Nw = 100
w_min = 0.0
w_max = 7.0
omega = np.linspace(w_min, w_max, Nw)

k = np.arange(lenk1)

intensity_sc1 = np.zeros((Nw, lenk1))
# intensity_sc2 = np.zeros((Nw, lenk1))
# intensity_sc3 = np.zeros((Nw, lenk1))

"""
when adding the magnetic field in experiments,
one domain was selected out
"""

for flag in range(lenk1):
    
    q1 = domain1[flag, 0]
    q2 = domain1[flag, 1]
    q3 = domain1[flag, 2]
    
    kx, ky, kz = cf.k12Tokxy(q1, q2, q3)
    intensity_sc1[:, flag] = GLSW.intensity(omega, kx, ky, kz)

# =============================================================================
# for flag in range(lenk1):
#     
#     q1 = domain2[flag, 0]
#     q2 = domain2[flag, 1]
#     q3 = domain2[flag, 2]
#     
#     kx, ky, kz = cf.k12Tokxy(q1, q2, q3)
#     intensity_sc2[:, flag] = GLSW.intensity(omega, kx, ky, kz)
#     
# for flag in range(lenk1):
#     
#     q1 = domain3[flag, 0]
#     q2 = domain3[flag, 1]
#     q3 = domain3[flag, 2]
#     
#     kx, ky, kz = cf.k12Tokxy(q1, q2, q3)
#     intensity_sc3[:, flag] = GLSW.intensity(omega, kx, ky, kz)
# =============================================================================
    
fname = 'cuts/K=' + str(K) + '/' + str(field) + 'T_intensity_glsw_d1.txt'
np.savetxt(fname, intensity_sc1)

# =============================================================================
# fname = 'intensity_2.txt'
# np.savetxt(fname, intensity_sc2)
# 
# fname = 'intensity_3.txt'
# np.savetxt(fname, intensity_sc3)
# =============================================================================
    
    
    
