#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 20:15:22 2019

@author: hao
"""

import numpy as np
import cofig as cf
import ex_info as exi

f0 = cf.f0
f1 = cf.f1
f2 = cf.f2
f3 = cf.f3
sub_idx = cf.sub_idx
delta_ij = cf.delta_ij
num_bond = cf.num_bond
num_sub = cf.num_sub

tHij = exi.tHij
thi = exi.thi


"""
obtain the bond depedent function in real space
there are four types, JJI, IIJ, JJJ, and III,
where JJI and IIJ are functions of (bond, flavor1(n), flavor2(n'), flavor3(m))
while JJJ and III are functions of (bond, flavor1(n), flavor2(n'))

"""

JJI = np.zeros((num_bond, 2, 2, 2), dtype=complex)
IIJ = np.zeros((num_bond, 2, 2, 2), dtype=complex)
JJJ = np.zeros((num_bond, 2, 2), dtype=complex)
III = np.zeros((num_bond, 2, 2), dtype=complex)

for bond in range(num_bond):

    THij = tHij[bond, :, :]

    for flavor1 in range(2):
        for flavor2 in range(2):

            f3nnp = f3[:, flavor1, flavor2]
            JJJ[bond, flavor1, flavor2] = f0 @ THij @ f3nnp
            III[bond, flavor1, flavor2] = f3nnp @ THij @ f0

            for flavor3 in range(2):
                f1m = f1[:, flavor3]
                f2nnp = f2[:, flavor1, flavor2]
                JJI[bond, flavor1, flavor2, flavor3] = f1m @ THij @ f2nnp
                IIJ[bond, flavor1, flavor2, flavor3] = f2nnp @ THij @ f1m
                

def phase_fun(q, vec, typ):
    """
    calculate the phase coming from Fourier transform
    typ 0: iii and jjj
    typ 1: iij
    typ 2: jji

    Parameters
    ----------
    q : np.array((3,))
        the momentum.
    vec : np.array((3,))
        the bond vector.
    typ : int
        type.

    Returns
    -------
    complex number
        the phase factor.

    """
    if typ == 0:                #iii or jjj
        return 1.0
    elif typ == 1:              #iij
        return np.exp(1j*q@vec)
    elif typ == 2:              #jji
        return np.exp(-1j*q@vec)
    else:
        print("wrong type number")
        return
    

def V_cubic_decay(q1, q2, q3, ubov1, ubov2, ubov3, ubovm1, ubovm2, ubovm3):
    
    """
    calculates the decay vertex functions: 
        V_decay = V2*beta^{\dagger}*beta^{\dagger}*beta

    Parameters
    ----------
    q1 : np.array((3,))
        input momentum.
    q2 : np.array((3,))
        intermediate momentum.
    q3 : np.array((3,))
        intermediate momentum.
    ubov{}: np.array((4*num_sub, 4*num_sub))
        Bogoliubov matrices for positive and negative momenta

    Returns
    -------
    V_decay : np.array((2*num_sub, 2*num_sub, 2*num_sub))  decay vertex
    """
    
    def Fcd_mat_symm(sub1, sub2, bond_vec, typ):
        
        Fc_mat_symm = np.zeros((2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                               dtype=complex)
        Fd_mat_symm = np.zeros((2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                               dtype=complex)
        
        
            
        
            
        
    

        
    
