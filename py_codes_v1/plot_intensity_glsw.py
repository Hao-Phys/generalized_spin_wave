#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:41:57 2019

@author: hao
"""

import numpy as np
#import matplotlib
#matplotlib.rcParams['text.usetex'] = True
from matplotlib import pyplot as plt

data1 = np.loadtxt('intensity_1.txt')
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
ax.set_ylabel('$\omega$ meV')
#plt.show()
plt.savefig('/home/hao/Desktop/int_test.png', dpi=100)