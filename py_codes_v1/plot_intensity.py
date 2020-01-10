#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:41:57 2019

@author: hao
"""

import numpy as np
#matplotlib.rcParams['text.usetex'] = True
from matplotlib import pyplot as plt
import sys

fname = sys.argv[1]
data1 = np.loadtxt(fname)

#data2 = np.loadtxt('intensity_2.txt')
#data3 = np.loadtxt('intensity_3.txt')


lenO = len(data1[:, 0])
lenk = len(data1[0, :])
intensity_sc = data1 

k = np.arange(lenk)
Nw = 100
w_min = 0.0
w_max = 7.0
omega = np.linspace(w_min, w_max, Nw)

K, O = np.meshgrid(k, omega)


fig, ax = plt.subplots()
cm = plt.pcolormesh(K, O, intensity_sc, cmap='jet', vmin=0, vmax=2.5, shading='gouraud')
plt.colorbar(cm)
plt.xticks([0, lenk/4-1, lenk/2-1, 3*lenk/4-1, lenk-1], \
           [-2.0, -1.0, -0.0, 1.0, 2.0]) 
plt.xlabel(r'$H$')
plt.ylabel(r'$\omega$ meV')
plt.show()
#plt.savefig('/home/hao/Desktop/int_test.png', dpi=100)
