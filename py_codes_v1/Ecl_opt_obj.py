#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 19:46:37 2019

@author: hao
"""

import numpy as np
import cofig as cf
import fun_generators as fg

# pass the lattice information
sub_idx  = cf.sub_idx
gamma_ij = cf.gamma_ij
num_bond = cf.num_bond
num_sub  = cf.num_sub
gamma_ijc = np.conj(gamma_ij)

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


def Energy_cl(v):
    
    funval = 0.0
    
    for bond in range(num_bond):
        id1 = sub_idx[bond, 0]
        id2 = sub_idx[bond, 1]
        
        alpha1 = v[4*id1]
        alpha2 = v[4*id1+1]
        theta  = v[4*id1+2]
        phi    = v[4*id1+3]
        
        alpha1_p = v[4*id2]
        alpha2_p = v[4*id2+1]
        theta_p  = v[4*id2+2]
        phi_p    = v[4*id2+3]
        
        S1p = fg.fun_sp(alpha1, alpha2, theta, phi)
        S1m = fg.fun_sm(alpha1, alpha2, theta, phi)
        S1z = fg.fun_sz(theta, phi)
        
        S2p = fg.fun_sp(alpha1_p, alpha2_p, theta_p, phi_p)
        S2m = fg.fun_sm(alpha1_p, alpha2_p, theta_p, phi_p)
        S2z = fg.fun_sz(theta_p, phi_p)
        
        funval += Jex[bond] * (0.5*(S1p*S2m+S1m*S2p) + Delta[bond]*S1z*S2z) \
                + Jpm[bond] * (gamma_ij[bond]*S1p*S2p + gamma_ijc[bond]*S1m*S2m) \
                - 0.5*1j*Jzpm[bond]*(gamma_ijc[bond]*S1p*S2z \
                                   - gamma_ij[bond]*S1m*S2z  \
                                   + gamma_ijc[bond]*S1z*S2p \
                                   - gamma_ij[bond]*S1z*S2m)
                                     
    for sub_lat in range(num_sub):
        theta = v[4*sub_lat+2]
        phi   = v[4*sub_lat+3]
        Sz = fg.fun_sz(theta, phi)
        Sz_sq = fg.fun_sz_sq(theta)        
        funval += - (D_ion*Sz_sq + h_ext*Sz)
        
    funval = np.real(funval)
    funval = float(funval)    
    return funval
        


                                     
     
        
        
        
        
        
    


