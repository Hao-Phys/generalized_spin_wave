#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 15:34:11 2019

@author: Hao
"""
import numpy as np
import cofig as cf
import GLSW
from matplotlib import pyplot as plt
#from matplotlib import cm

K = 0.1
inputpath1 = 'cuts/K=' + str(K) + '/path_domain1.dat'
domain1 = np.loadtxt(inputpath1)

lenk1 = len(domain1[:, 0])

# the energy range
Nw = 100
w_min = 0.0
w_max = 7.0
omega = np.linspace(w_min, w_max, Nw)

k = np.arange(lenk1)

intensity_sc1 = np.zeros((Nw, lenk1))

for flag in range(lenk1):
    
    q1 = domain1[flag, 0]
    q2 = domain1[flag, 1]
    q3 = domain1[flag, 2]
    
    kx, ky, kz = cf.k12Tokxy(q1, q2, q3)
    intensity_sc1[:, flag] = GLSW.intensity(omega, kx, ky, kz)
    
K, O = np.meshgrid(k, omega)
c = plt.pcolormesh(K, O, intensity_sc1, cmap='viridis')
plt.show()


    
    
    