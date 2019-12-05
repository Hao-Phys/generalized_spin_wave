#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 21:24:42 2019

@author: hao
"""


import numpy as np
import sys
"""
model parameters
"""
# accept input files
# run as:  python xxx.py input.txt

inFile = sys.argv[1]   
paras = np.loadtxt(inFile)

#paras = np.loadtxt('input1.txt')

kelvin_to_meV = 11.602

# model parameters in meV

J1      = paras[0]/kelvin_to_meV
Delta1  = paras[1]
J2      = paras[2]/kelvin_to_meV
Delta2  = paras[3]
J3      = paras[4]/kelvin_to_meV
Delta3  = paras[5]
D_ion   = paras[6]/kelvin_to_meV  
J1pm    = paras[7]/kelvin_to_meV/2
J3pm    = paras[8]/kelvin_to_meV
J1zpm   = paras[9]/kelvin_to_meV 
J3zpm   = paras[10]/kelvin_to_meV
J0p     = paras[11]/kelvin_to_meV
Delta0p = paras[12]
J1p     = paras[13]/kelvin_to_meV 
Delta1p = paras[14]
J2ap    = paras[15]/kelvin_to_meV
Delta2ap= paras[16]
J2bp    = paras[17]/kelvin_to_meV
Delta2bp= paras[18]
field   = paras[19]

print("field is:  ", field)
# directory for each field for later manipulation
path = 'gs_info/h=' + str(field) +'T/'


# the external field
g_factor = 3.2
mu_B = 5.788381E-2
h_ext = g_factor*mu_B*field

# the lattice information

num_sub  = 4
num_bond = 88

# broadening factor

broadening = 0.12

# the A-matrix

A_mat = np.zeros((4*num_sub, 4*num_sub))

for flag in range(2*num_sub):
    A_mat[flag, flag] = 1
    
for flag in range(2*num_sub, 4*num_sub):
    A_mat[flag, flag] = -1
    
### the bond information ###
 
a1 = np.array([[1, 0, 0]])
a2 = np.array([[0, 1, 0]])
a3_true = np.array([[0, 0, 1]])
a3 = a3_true - (a2-a1/4)
 
sub_idx  = np.zeros((num_bond, 2))
delta_ij = np.zeros((3, num_bond))
gamma_ij = np.zeros(num_bond, dtype = complex)


for flag in range(4):    
    # intra-layer bonds
    
    # nearest neighbor bonds
    sub_idx[flag*3, 0]  = flag
    sub_idx[flag*3, 1]  = np.mod(flag+1, 4)
    delta_ij[:, flag*3] = 0.25*a1 - a2
    gamma_ij[flag*3]    = 1.0
    
    sub_idx[flag*3+1, 0]  = flag
    sub_idx[flag*3+1, 1]  = np.mod(flag+1, 4)
    delta_ij[:, flag*3+1] = 0.25*a1
    gamma_ij[flag*3+1]    = np.exp(-1j*2*np.pi/3)
    
    sub_idx[flag*3+2, 0]  = flag
    sub_idx[flag*3+2, 1]  = flag
    delta_ij[:, flag*3+2] = a2
    gamma_ij[flag*3+2]    = np.exp(1j*2*np.pi/3)
    
    # second nearest neighbor bonds
    sub_idx[flag*3+12, 0]  = flag
    sub_idx[flag*3+12, 1]  = np.mod(flag+2, 4)
    delta_ij[:, flag*3+12] = 0.5*a1 - a2 
    
    sub_idx[flag*3+13, 0]  = flag
    sub_idx[flag*3+13, 1]  = np.mod(flag+1, 4)
    delta_ij[:, flag*3+13] = 0.25*a1 + a2
    
    sub_idx[flag*3+14, 0]  = flag
    sub_idx[flag*3+14, 1]  = np.mod(flag-1, 4)
    delta_ij[:, flag*3+14] = -0.25*a1 + 2*a2
    
    # third nearest neighbor bonds
    sub_idx[flag*3+24, 0]  = flag
    sub_idx[flag*3+24, 1]  = np.mod(flag+2, 4)
    delta_ij[:, flag*3+24] = 0.5*a1 - 2*a2
    gamma_ij[flag*3+24]    = 1.0
    
    sub_idx[flag*3+25, 0]  = flag
    sub_idx[flag*3+25, 1]  = np.mod(flag+2, 4)
    delta_ij[:, flag*3+25] = 0.5*a1
    gamma_ij[flag*3+25]    = np.exp(-1j*2*np.pi/3)
    
    sub_idx[flag*3+26, 0]  = flag
    sub_idx[flag*3+26, 1]  = flag
    delta_ij[:, flag*3+26] = 2.0*a2
    gamma_ij[flag*3+26]    = np.exp(1j*2*np.pi/3)
    
    # inter-layer bonds
    
    # J0' bonds
    sub_idx[36+flag, 0]    = flag
    sub_idx[36+flag, 1]    = np.mod(flag+1, 4)
    delta_ij[:, 36+flag]   = a3
    
    # J1' bonds
    # bond 1
    sub_idx[40+flag*6, 0]  = flag
    sub_idx[40+flag*6, 1]  = np.mod(flag+2, 4)
    delta_ij[:, 40+flag*6] = a1/4 - a2 + a3
    
    # bond 2
    sub_idx[41+flag*6, 0]  = flag
    sub_idx[41+flag*6, 1]  = np.mod(flag+2, 4)
    delta_ij[:, 41+flag*6] = a1/4 + a3
    
    # bond 3
    sub_idx[42+flag*6, 0]  = flag
    sub_idx[42+flag*6, 1]  = np.mod(flag+1, 4)
    delta_ij[:, 42+flag*6] = a2 + a3
    
    # bond 4
    sub_idx[43+flag*6, 0]  = flag
    sub_idx[43+flag*6, 1]  = flag
    delta_ij[:, 43+flag*6] = a3 - (a1/4 - a2)
    
    # bond 5
    sub_idx[44+flag*6, 0]  = flag
    sub_idx[44+flag*6, 1]  = flag
    delta_ij[:, 44+flag*6] = a3 - a1/4 
    
    # bond 6
    sub_idx[45+flag*6, 0]  = flag
    sub_idx[45+flag*6, 1]  = np.mod(flag+1, 4)
    delta_ij[:, 45+flag*6] = a3 - a2
    
    # J2a' bonds
    # bond 1
    sub_idx[64+flag*3, 0]  = flag
    sub_idx[64+flag*3, 1]  = np.mod(flag+2, 4)
    delta_ij[:, 64+flag*3] = a3 + a1/4 - 2*a2
    
    # bond 2 
    sub_idx[65+flag*3, 0]  = flag
    sub_idx[65+flag*3, 1]  = np.mod(flag+2, 4)
    delta_ij[:, 65+flag*3] = a3 + a1/4 + a2
    
    # bond 3
    sub_idx[66+flag*3, 0]  = flag
    sub_idx[66+flag*3, 1]  = np.mod(flag+3, 4)
    delta_ij[:, 66+flag*3] = a3 + a2 - a1/2
    
    # J2b' bonds
    # bond 1
    sub_idx[76+flag*3, 0]  = flag
    sub_idx[76+flag*3, 1]  = np.mod(flag+3, 4)
    delta_ij[:, 76+flag*3] = a3 + a1/2 - a2
    
    # bond 2
    sub_idx[77+flag*3, 0]  = flag
    sub_idx[77+flag*3, 1]  = flag
    delta_ij[:, 77+flag*3] = a3 - (a1/4 - 2*a2)
    
    # bond 3
    sub_idx[78+flag*3, 0]  = flag
    sub_idx[78+flag*3, 1]  = flag
    delta_ij[:, 78+flag*3] = a3 - (a1/4 + a2)
    
sub_idx = sub_idx.astype(int)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
    
    


