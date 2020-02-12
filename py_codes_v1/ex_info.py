#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : ex_info.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.06.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:14:49 2019

@author: hao
"""

import numpy as np
import cofig as cf
import SU3rot as rot

# pass the lattice information
sub_idx  = cf.sub_idx
gamma_ij = cf.gamma_ij
num_bond = cf.num_bond
num_sub  = cf.num_sub
#gamma_ijc = np.conj(gamma_ij)

# the exchange tensor in the local reference frame
tHij = np.zeros((num_bond, 8, 8))
# on-site term in the local reference frame
thi = np.zeros((num_sub, 8))

Jex = np.zeros((num_bond, 1))
Jpm = np.zeros((num_bond, 1))
Jzpm = np.zeros((num_bond, 1))
Delta = np.zeros((num_bond, 1))

# pass the parameters
Jex[0:12] = cf.J1
Jex[12:24] = cf.J2
Jex[24:36] = cf.J3
Jex[36:40] = cf.J0p
Jex[40:64] = cf.J1p
Jex[64:76] = cf.J2ap
Jex[76:88] = cf.J2bp
Jpm[0:12] = cf.J1pm
Jpm[12:24] = 0.0
Jpm[24:36] = cf.J3pm
Jpm[36:88] = 0.0
Jzpm[0:12] = cf.J1zpm
Jzpm[12:24] = 0.0
Jzpm[24:36] = cf.J3zpm
Jzpm[36:88] = 0.0
Delta[0:12] = cf.Delta1
Delta[12:24] = cf.Delta2
Delta[24:36] = cf.Delta3
Delta[36:40] = cf.Delta0p
Delta[40:64] = cf.Delta1p
Delta[64:76] = cf.Delta2ap
Delta[76:88] = cf.Delta2bp
D_ion = cf.D_ion
h_ext = cf.h_ext

for bond in range(num_bond):
    # exchange tensor in the global frame
    tmp = np.zeros((8, 8))
    tmp[0, 0] = Jex[bond] + 2*np.real(gamma_ij[bond])*Jpm[bond]
    tmp[0, 1] = -2*np.imag(gamma_ij[bond])*Jpm[bond]
    tmp[0, 2] = -np.imag(gamma_ij[bond])*Jzpm[bond]
    tmp[1, 0] = -2*np.imag(gamma_ij[bond])*Jpm[bond]
    tmp[1, 1] = Jex[bond] - 2*np.real(gamma_ij[bond])*Jpm[bond]
    tmp[1, 2] = np.real(gamma_ij[bond])*Jzpm[bond] 
    tmp[2, 0] = -np.imag(gamma_ij[bond])*Jzpm[bond]
    tmp[2, 1] = np.real(gamma_ij[bond])*Jzpm[bond]  
    tmp[2, 2] = Jex[bond]*Delta[bond]

    # if bond in range(25, 28):
        # tmpp = tmp[:3, :3]
        # trans = np.array([[0.5, 0.5, 0],
                          # [0.5/(1j), -(0.5/1j), 0],
                          # [0, 0, 1.0]])
        # tmp3 = trans.conj().T @ tmpp @ trans
        # print(tmp3)
    sub1 = sub_idx[bond, 0]
    sub2 = sub_idx[bond, 1]
    
    R1 = rot.R_mat[sub1, :, :]
    R2 = rot.R_mat[sub2, :, :]
    tHij[bond, :, :] = R1 @ tmp @ R2.T
    
for sublat in range(num_sub):
    tmp1 = np.zeros((8, 1))
    tmp1[2] = -h_ext
    tmp1[7] = -1.0/np.sqrt(3)*D_ion
    
    RR = rot.R_mat[sublat, :, :]
    thi[sublat, :] = (RR @ tmp1).reshape(8)
    
    
    
