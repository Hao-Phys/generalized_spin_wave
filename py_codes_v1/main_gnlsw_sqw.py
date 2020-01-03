#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 21:35:51 2020

@author: hao
"""
"""
To run this program
   python3 main_gnlsw_sqw.py <input.txt> <selfE.txt>
   where <input.txt> is the input file, and 
   <selfE.txt> is the file stores the self energies (on-shell)

"""
import numpy as np
import cofig as cf
import sys
import GNLSW_selfE

field = cf.field

inFile = sys.argv[2]
dirname = inFile.rstrip(inFile.rsplit('/')[-1])


data = np.loadtxt(inFile)
lenf = len(data[:, 0])

num_sub = cf.num_sub

# the energy range
Nw = 100
w_min = 0.0
w_max = 7.0
omega = np.linspace(w_min, w_max, Nw)

intensity_nlsw = np.zeros((Nw, lenf))

for flag in range(lenf):
    q1 = data[flag, 1]
    q2 = data[flag, 2]
    q3 = data[flag, 3]
    q = np.array([q1, q2, q3])
    selfE_re = np.zeros([2*num_sub])
    selfE_im = np.zeros([2*num_sub])
    
    for band in range(2*num_sub):
        selfE_re[band] = data[flag, 2*band+4] + data[flag, band+20]
        selfE_im[band] = data[flag, 2*band+5]
    
    kx, ky, kz = cf.k12Tokxy(q1, q2, q3)
    intensity_nlsw[:, flag] = GNLSW_selfE.intensity(omega, kx, ky, kz, \
                                                    selfE_re, selfE_im)
    
fname = dirname + str(field) + 'T_intensity_nlsw.txt'
np.savetxt(fname, intensity_nlsw) 
    
     
        
