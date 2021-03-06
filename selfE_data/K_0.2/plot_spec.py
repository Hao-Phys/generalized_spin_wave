#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 21:25:38 2020

@author: hao
"""

import numpy as np
# =============================================================================
# import cofig as cf
# import GLSW
# =============================================================================
import sys
from matplotlib import pyplot as plt

inFile = sys.argv[1]
#dirname = inFile.rstrip(inFile.rsplit('/')[-1])
data = np.loadtxt(inFile)
lenf = len(data[:, 0])
k = np.arange(lenf)

# =============================================================================
# data = np.loadtxt(inFile)
# lenf = len(data[:, 0])
# k = np.arange(lenf)
# 
# num_sub = cf.num_sub
# 
# ek = np.zeros((lenf, 2*num_sub))
# omega = np.zeros((lenf, 2*num_sub))
# Gamma = np.zeros((lenf, 2*num_sub))
# 
# for flag in range(lenf):
#     q1 = data[flag, 1]
#     q2 = data[flag, 2]
#     q3 = data[flag, 3]
#     q = np.array([q1, q2, q3])
#     energy, tmp = GLSW.eigensystem(q)
#     ek[flag, :] = energy[:2*num_sub] 
# 
#     for band in range(2*num_sub):
#         spec_nlsw = energy[band] + data[flag, 2*band+4] + data[flag, band+20]
#         omega[flag, band] = spec_nlsw
#         decay_rate = data[flag, 2*band+5]
#         Gamma[flag, band] = decay_rate
# 
# 
# ekmax = omega + Gamma
# ekmin = omega - Gamma
# 
# fig, ax = plt.subplots()
# =============================================================================

plt.plot(k, data[:, 11])
    
# =============================================================================
# plt.xlim(0, lenf-1)
# 
# =============================================================================
#plt.ylim(-0.3, 0)
plt.show()




