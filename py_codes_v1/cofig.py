#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : cofig.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.27.2020
# Last Modified Date: 01.27.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
import numpy as np
import sys
import os

"""
model parameters
"""

# accept input fiels
# run the program as python{} xxx.py input.txt

inFile = sys.argv[1]
#inFile = 'input.txt'
paras = np.loadtxt(inFile)

kelvin_to_meV = 11.602

# model parameters in meV

J1      = paras[0]/kelvin_to_meV
Delta1  = paras[1]
J2      = paras[2]/kelvin_to_meV;
Delta2  = paras[3]
J3      = paras[4]/kelvin_to_meV
Delta3  = paras[5]
D_ion   = paras[6]/kelvin_to_meV  
J1pm    = paras[7]/kelvin_to_meV/2.0
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

print("the external field is: ", field)
# directory for each field for later manipulation
cwd = os.getcwd()
path = cwd + '/gs_info/h=' + str(field) + 'T/'
# the external field

g_factor = 3.2
mu_B = 5.788381E-2
h_ext = g_factor*mu_B*field

# the lattice information

num_sub  = 4
num_bond = 88

# broadening factor

broadening = 0.06

# convergence factor, can be adjusted to the specific problem
convergence = 0.1

# the A-matrix

A_mat = np.zeros((4*num_sub, 4*num_sub))

for flag in range(2*num_sub):
    A_mat[flag, flag] = 1
    
for flag in range(2*num_sub, 4*num_sub):
    A_mat[flag, flag] = -1
    
### the bond information ###
 
a1 = np.array([[1.0, 0, 0]])
a2 = np.array([[0, 1.0, 0]])
a3_true = np.array([[0, 0, 1.0]])
a3 = a3_true - (a2-a1/4)
 
sub_idx  = np.zeros((num_bond, 2))
delta_ij = np.zeros((3, num_bond))
gamma_ij = np.zeros((num_bond, 1), dtype = complex)


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
    

# operators information
f0 = np.zeros(8)
f1 = np.zeros((8, 2), dtype = complex)
f2 = np.zeros((8, 2, 2), dtype = complex)
f3 = np.zeros((8, 2, 2), dtype = complex)

f0[2] = 1.0
f0[7] = 1.0/np.sqrt(3.0)

# 0 for the -1 boson, 1 for the 0 boson
f1[0, 1] = 1.0/np.sqrt(2.0)
f1[1, 1] = -1j/np.sqrt(2.0)
f1[3, 1] = 1j/np.sqrt(2.0)
f1[4, 1] = -1.0/np.sqrt(2.0)
f1[5, 0] = 1.0
f1[6, 0] = -1j

f2[0, 0, 1] = 0.5/np.sqrt(2.0)
f2[0, 1, 0] = 0.5/np.sqrt(2.0)
f2[1, 0, 1] = 0.5*1j/np.sqrt(2.0)
f2[1, 1, 0] = -0.5*1j/np.sqrt(2.0)
f2[2, 0, 0] = -1.0
f2[2, 1, 1] = -0.5
f2[3, 0, 1] = 0.5*1j/np.sqrt(2.0)
f2[3, 1, 0] = -0.5*1j/np.sqrt(2.0)
f2[4, 0, 1] = 0.5/np.sqrt(2.0)
f2[4, 1, 0] = 0.5/np.sqrt(2.0)
f2[7, 1, 1] = -1.5/np.sqrt(3.0)

f3[0, 0, 1] = -1.0/(2.0*np.sqrt(2.0))
f3[0, 1, 1] = -1.0/(2.0*np.sqrt(2.0))
f3[1, 0, 1] = 1j/(2.0*np.sqrt(2.0))
f3[1, 1, 1] = 1j/(2.0*np.sqrt(2.0))
f3[3, 0, 1] = -1j/(2.0*np.sqrt(2.0))
f3[3, 1, 1] = -1j/(2.0*np.sqrt(2.0))
f3[4, 0, 1] = 1.0/(2.0*np.sqrt(2.0))
f3[4, 1, 1] = 1.0/(2.0*np.sqrt(2.0))
f3[5, 0, 0] = -0.5
f3[5, 1, 0] = -0.5
f3[6, 0, 0] = 0.5*1j
f3[6, 1, 0] = 0.5*1j


# change of basis functions
def kxyTok12(kx, ky, kz):
    k1 = (kx + np.sqrt(3.0)*ky)/np.pi
    k2 = (np.sqrt(3.0)*ky - kx)/(4.0*np.pi)
    k3 = 2.0*kz/np.pi
    
    return k1, k2, k3

def k12Tokxy(k1, k2, k3):
    kx = np.pi * (0.5*k1 - 2.0*k2)
    ky = np.pi/np.sqrt(3) * (0.5*k1 + 2*k2)
    kz = 2.0*np.pi*k3/4.0
    
    return kx, ky, kz
# form factor

def formfactor(q):
    """
    the form factor of Fe2+ inons
    """
    
    A  = 0.0263
    aa = 34.9597
    B  = 0.3668
    bb = 15.9435 
    C  = 0.6188	
    cc = 5.5935
    D  = -0.0119
    
    s_sq = ((q[0]/(4.05))**2 + (q[1]/(4.05))**2 \
    + (q[2]/6.76)**2)/(4.0*np.pi)**2;
            
    ff= (A*np.exp(-aa*s_sq) + B*np.exp(-bb*s_sq) + C*np.exp(-cc*s_sq)+ D)**2
    
    return ff

def projector(qx, qy, qz, sqw_mat):
    """
    project out parallel components
    """
    #len_omega = sqw_mat.shape[0]
    #inten = np.zeros(len_omega)
    inten = 0.0
    q = np.array([qx, qy, qz])
    norm_sq = qx**2 + qy**2 + qz**2
    
    if norm_sq == 0:
        inten = sqw_mat[:, 0, 0] + sqw_mat[:, 1, 1] + sqw_mat[:, 2, 2]
    else:
        
# =============================================================================
        """
        projector using the numpy ufunc, currently does not win over the loop
        below.
        
        """
#
#       project_mat = np.array([[1.0-qx**2/norm_sq, -qx*qy/norm_sq, -qx*qz/norm_sq], \
#                                 [-qy*qx/norm_sq, 1.0-qy**2/norm_sq, -qy*qz/norm_sq], \
#                                 [-qz*qx/norm_sq, -qz*qy/norm_sq, 1.0-qz**2/norm_sq]])
#     
#         inten_mat = project_mat[None, :, :] * sqw_mat
#         inten = np.sum(np.sum(inten_mat, axis=2), axis=1)
#         
# =============================================================================
        for mu0 in range(3):
            for nu0 in range(3):
            
                if (mu0 == nu0):
                   inten += (1.0 - q[mu0]*q[nu0]/norm_sq) * sqw_mat[:, mu0, nu0]
                else:
                   inten += - q[mu0]*q[nu0]/norm_sq * sqw_mat[:, mu0, nu0]
                
    return inten
    
    
    



    
    
    
