#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:24:29 2019

@author: Hao
"""

import numpy as np
import cofig as cf
import ex_info as exi
import GLSW

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

def Fa_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):
    """
    The Fa function defined in the note
    """
    
    U21_q1 = Ubov1[2*num_sub:, :2*num_sub]
    U11_q2 = Ubov2[:2*num_sub, :2*num_sub]
    U11_q3 = Ubov3[:2*num_sub, :2*num_sub]
    phase = phase_fun(q3, vec, typ)
    value = U21_q1[X1, band1] * U11_q2[X2, band2] * U11_q3[X3, band3] * phase
    
    return value

def Fa_symm(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):
    
    """
    symmetrized Fa function: permutation of q1, q2, q3
    the band index and the order of Bogliubov coefficients
    change accordingly
    
    """
    value = 0.0
    
    tmp = Fa_fun(X1, X2, X3, \
                 band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
                 vec, typ)
    value += tmp
    
    tmp = Fa_fun(X1, X2, X3, \
                 band1, band3, band2, Ubov1, Ubov3, Ubov2, q1, q3, q2, \
                 vec, typ)
    value += tmp
        
    tmp = Fa_fun(X1, X2, X3, \
                 band2, band1, band3, Ubov2, Ubov1, Ubov3, q2, q1, q3, \
                 vec, typ)
    value += tmp
    
    tmp = Fa_fun(X1, X2, X3, \
                 band2, band3, band1, Ubov2, Ubov3, Ubov1, q2, q3, q1, \
                 vec, typ)
    value += tmp
    
    tmp = Fa_fun(X1, X2, X3, \
                 band3, band2, band1, Ubov3, Ubov2, Ubov1, q3, q2, q1, \
                 vec, typ)
    value += tmp
    
    tmp = Fa_fun(X1, X2, X3, \
                 band3, band1, band2, Ubov3, Ubov1, Ubov2, q3, q1, q2, \
                 vec, typ)
    value += tmp
    
    return value # avoid overcouting
    
    
def Fb_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):

    U11_mq1 = Ubov1[:2*num_sub, :2*num_sub].conj()
    U21_mq2 = Ubov2[2*num_sub:, :2*num_sub].conj()
    U21_mq3 = Ubov3[2*num_sub:, :2*num_sub].conj()
    phase = phase_fun(-q3, vec, typ)
    value = U11_mq1[X1, band1] * U21_mq2[X2, band2] * U21_mq3[X3, band3] * phase
    
    return value

def Fb_symm(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):
    
    value = 0.0
    
    tmp = Fb_fun(X1, X2, X3, \
                 band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
                 vec, typ)
    value += tmp
    
    tmp = Fb_fun(X1, X2, X3, \
                 band1, band3, band2, Ubov1, Ubov3, Ubov2, q1, q3, q2, \
                 vec, typ)
    value += tmp
        
    tmp = Fb_fun(X1, X2, X3, \
                 band2, band1, band3, Ubov2, Ubov1, Ubov3, q2, q1, q3, \
                 vec, typ)
    value += tmp
    
    tmp = Fb_fun(X1, X2, X3, \
                 band2, band3, band1, Ubov2, Ubov3, Ubov1, q2, q3, q1, \
                 vec, typ)
    value += tmp
    
    tmp = Fb_fun(X1, X2, X3, \
                 band3, band2, band1, Ubov3, Ubov2, Ubov1, q3, q2, q1, \
                 vec, typ)
    value += tmp
    
    tmp = Fb_fun(X1, X2, X3, \
                 band3, band1, band2, Ubov3, Ubov1, Ubov2, q3, q1, q2, \
                 vec, typ)
    value += tmp
    
    return value  # avoid overcouting
    
    
def Fc_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):

    U11_mq1 = Ubov1[:2*num_sub, :2*num_sub].conj()
    U21_mq1 = Ubov1[2*num_sub:, :2*num_sub].conj()
    U21_mq2 = Ubov2[2*num_sub:, :2*num_sub].conj()
    U11_q3 = Ubov3[:2*num_sub, :2*num_sub]
    U21_q3 = Ubov3[2*num_sub:, :2*num_sub]
    phase1 = phase_fun(q3, vec, typ)
    phase2 = phase_fun(-q1, vec, typ)
    phase3 = phase_fun(-q2, vec, typ)
    value = U11_mq1[X1, band1]*U21_mq2[X2, band2]*U11_q3[X3, band3]*phase1 \
          + U21_q3[X1, band3]*U21_mq2[X2, band2]*U21_mq1[X3, band1]*phase2 \
          + U11_mq1[X1, band1]*U11_q3[X2, band3]*U21_mq2[X3, band2]*phase3
          
    return value

def Fc_symm(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):
    
    value = 0.0
    
    tmp = Fc_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ)
    value += tmp
    
    tmp = Fc_fun(X1, X2, X3, \
           band2, band1, band3, Ubov2, Ubov1, Ubov3, q2, q1, q3, \
           vec, typ)
    value += tmp
    
    return value # avoid overcouting



def Fd_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):

    U11_mq1 = Ubov1[:2*num_sub, :2*num_sub].conj()
    U21_mq1 = Ubov1[2*num_sub:, :2*num_sub].conj()
    U11_q2 = Ubov2[:2*num_sub, :2*num_sub]
    U21_q2 = Ubov2[2*num_sub:, :2*num_sub]
    U11_q3 = Ubov3[:2*num_sub, :2*num_sub]
    U21_q3 = Ubov3[2*num_sub:, :2*num_sub]
    phase1 = phase_fun(q3, vec, typ)
    phase2 = phase_fun(-q1, vec, typ)
    value = U11_mq1[X1, band1]*U11_q2[X2, band2]*U11_q3[X3, band3]*phase1 \
          + U21_q3[X1, band3]*U11_q2[X2, band2]*U21_mq1[X3, band1]*phase2 \
          + U21_q2[X1, band2]*U21_mq1[X2, band1]*U11_q3[X3, band3]*phase1
          
    return value


def Fd_symm(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ):
    
    value = 0.0
    
    tmp = Fd_fun(X1, X2, X3, \
           band1, band2, band3, Ubov1, Ubov2, Ubov3, q1, q2, q3, \
           vec, typ)
    value += tmp
    
    tmp = Fd_fun(X1, X2, X3, \
           band1, band3, band2, Ubov1, Ubov3, Ubov2, q1, q3, q2, \
           vec, typ)
    
    value += tmp
    
    return value

# the vertex functions
def V1_cubic(band1, band2, band3, q1, q2, q3):
    
    V1 = 0.0
    
    tmp, Ubov1 = GLSW.eigensystem(q1)
    tmp, Ubov2 = GLSW.eigensystem(q2)
    tmp, Ubov3 = GLSW.eigensystem(q3)
    
    tmp, Ubovm1 = GLSW.eigensystem(-q1)
    tmp, Ubovm2 = GLSW.eigensystem(-q2)
    tmp, Ubovm3 = GLSW.eigensystem(-q3)

    for N in range(2):
      for Np in range(2):
        
 
        for bond in range(1):
            
        
           bond_vec = delta_ij[:, bond]
           sub1 = sub_idx[bond, 0]
           sub2 = sub_idx[bond, 1]
           
           # N, Np, and M are flavor indices defined in the note
           In = 4*N + sub1
           Jn = 4*N + sub2
           Inp = 4*Np + sub1
           Jnp = 4*Np + sub2 
                       
           tFa1 = Fa_symm(Jn, Jn, Jnp, \
                          band1, band2, band3, Ubov1, Ubov2, Ubov3, \
                          q1, q2, q3, bond_vec, 0)
            
           tFb1 = Fb_symm(Jn, Jn, Jnp, \
                          band1, band2, band3, Ubovm1, Ubovm2, Ubovm3, \
                          -q1, -q2, -q3, bond_vec, 0)
            
           tFa2 = Fa_symm(In, In, Inp, \
                          band1, band2, band3, Ubov1, Ubov2, Ubov3, \
                          q1, q2, q3, bond_vec, 0)
            
           tFb2 = Fb_symm(In, In, Inp, \
                          band1, band2, band3, Ubovm1, Ubovm2, Ubovm3, \
                          -q1, -q2, -q3, bond_vec, 0)
            
           V1 += JJJ[bond, N, Np]*tFa1 + (JJJ[bond, N, Np]*tFb1).conj() \
                +III[bond, N, Np]*tFa2 + (III[bond, N, Np]*tFb2).conj()
                 
           for M in range(2):
                
                Im = 4*M + sub1
                Jm = 4*M + sub2
                
                tFa1 = Fa_symm(Jn, Jnp, Im, \
                               band1, band2, band3, Ubov1, Ubov2, Ubov3, \
                               q1, q2, q3, bond_vec, 2)
                
                tFb1 = Fb_symm(Jn, Jnp, Im, \
                               band1, band2, band3, Ubovm1, Ubovm2, Ubovm3, \
                               -q1, -q2, -q3, bond_vec, 2)
                
                tFa2 = Fa_symm(In, Inp, Jm, \
                               band1, band2, band3, Ubov1, Ubov2, Ubov3, \
                               q1, q2, q3, bond_vec, 1)
                
                tFb2 = Fb_symm(In, Inp, Jm, \
                               band1, band2, band3, Ubovm1, Ubovm2, Ubovm3, \
                               -q1, -q2, -q3, bond_vec, 1)
                
                V1 += 2.0*(JJI[bond, N, Np, M]*tFa1 \
                           + (JJI[bond, N, Np, M]*tFb1).conj() \
                         + IIJ[bond, N, Np, M]*tFa2 \
                           + (IIJ[bond, N, Np, M]*tFb2).conj())
                    
        for sublat in range(num_sub):
            
            In = 4*N + sublat
            Inp = 4*Np + sublat
            factor1 = thi[sublat] @ f3[:, N, Np]
            print(thi[sublat])
           
            factor2 = thi[sublat] @ f3[:, N, Np].conj()
            tFa = Fa_symm(In, In, Inp, \
                          band1, band2, band3, Ubov1, Ubov2, Ubov3, \
                          q1, q2, q3, 0.0, 0)
            tFb = Fb_symm(In, In, Inp, \
                          band1, band2, band3, Ubovm1, Ubovm2, Ubovm3, \
                          -q1, -q2, -q3, 0.0, 0)
            
            V1 += factor1*tFa + factor2*tFb.conj()
            
    return V1
            
        
                    
                
                
                
                


 
