#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:24:29 2019

@author: Hao
"""

import numpy as np
import cofig as cf
import ex_info as exi

f0 = cf.f0
f1 = cf.f1
f2 = cf.f2
f3 = cf.f3
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
                JJI[bond, flavor1, flavor2, flavor3] = f1m @ f2nnp
                IIJ[bond, flavor1, flavor2, flavor3] = f2nnp @ f1m
                

# definition of the phase function
def phase_fun(q, vec, typ):
    if typ == 0:                #iii or jjj
        return 1.0
    elif typ == 1:              #iij
        return np.exp(1j*q@vec) 
    elif typ == 2:              #jji
        return np.exp(-1j*q@vec)
    else:
        print("wrong type number")
        return

                
# the basic ingredients
                                
def Fa_fun(band1, band2, band3, X1, X2, X3, \
           U21_q1, U11_q2, U11_q3, q3, vec, typ):
        
    phase = phase_fun(q3, vec, typ)
    value = U21_q1[X1, band1] * U11_q2[X2, band2] * U11_q3[X3, band3] * phase
    return value
    
def Fb_fun(band1, band2, band3, X1, X2, X3, \
           U11_mq1, U21_mq2, U21_mq3, q3, vec, typ):
    
    phase = phase_fun(-q3, vec, typ)
    value = U11_mq1[X1, band1] * U21_mq2[X2, band2] * U21_mq3[X3, band3] * phase
    return value


    
        
            
            
            






