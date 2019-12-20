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
III_onsite = np.zeros((num_sub, 2, 2), dtype=complex)
III_onsite_c = np.zeros((num_sub, 2, 2), dtype=complex)


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
 
               
for sublat in range(num_sub):
    
    for flavor1 in range(2):
        for flavor2 in range(2):
            
            III_onsite[sublat, flavor1, flavor2] = thi[sublat, :] \
                @ f3[:, flavor1, flavor2]
            III_onsite_c[sublat, flavor1, flavor2] = thi[sublat, :] \
                @ f3[:, flavor1, flavor2].conj()
                

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
    
#    u11_1 = ubov1[:2*num_sub, :2*num_sub]
#    u21_1 = ubov1[2*num_sub:, :2*num_sub]
#    u11_2 = ubov2[:2*num_sub, :2*num_sub]
#    u21_2 = ubov2[2*num_sub:, :2*num_sub]
    u11_3 = ubov3[:2*num_sub, :2*num_sub]
    u21_3 = ubov3[2*num_sub:, :2*num_sub]
    
    u11_m1 = ubovm1[:2*num_sub, :2*num_sub]
    u21_m1 = ubovm1[2*num_sub:, :2*num_sub]
    u11_m2 = ubovm2[:2*num_sub, :2*num_sub]
    u21_m2 = ubovm2[2*num_sub:, :2*num_sub]
#    u11_m3 = ubovm3[:2*num_sub, :2*num_sub]
#    u21_m3 = ubovm3[2*num_sub:, :2*num_sub]
    
    def Fcd_mat_symm(sub1, sub2, bond_vec, typ):
        
        
        phase_1 = phase_fun(q1, bond_vec, typ)
        phase_2 = phase_fun(q2, bond_vec, typ)
        phase_3 = phase_fun(q3, bond_vec, typ)
        phase_m1 = phase_fun(-q1, bond_vec, typ)
        phase_m2 = phase_fun(-q2, bond_vec, typ)
        phase_m3 = phase_fun(-q3, bond_vec, typ)
        
        if (typ == 0):
            
           Fc_mat_symm = np.zeros((2, 2, 2*num_sub, 2*num_sub, \
                                           2*num_sub), dtype=complex)
           Fd_mat_symm = np.zeros((2, 2, 2*num_sub, 2*num_sub, \
                                           2*num_sub), dtype=complex)
           for N in range(2):
               X1 = num_sub*N + sub1
               X2 = X1
               
               for Np in range(2):
                   
                   X3 = num_sub*Np + sub1
                   
                   tmp1 = np.outer(u11_m1[X1, :], u21_m2[X2, :])[:, :, None] \
                           * u11_3[X3, :][None, None, :] * phase_3 \
                           + np.outer(u21_m1[X3, :], u21_m2[X2, :])[:, :, None] \
                           * u21_3[X1, :][None, None, :] * phase_m1 \
                           + np.outer(u11_m1[X1, :], u21_m2[X3, :])[:, :, None] \
                           * u11_3[X2, :][None, None, :] * phase_m2 \
                           + np.outer(u21_m1[X2, :], u11_m2[X1, :])[:, :, None] \
                           * u11_3[X3, :][None, None, :] * phase_3 \
                           + np.outer(u21_m1[X2, :], u21_m2[X3, :])[:, :, None] \
                           * u21_3[X1, :][None, None, :] * phase_m2 \
                           + np.outer(u21_m1[X3, :], u11_m2[X1, :])[:, :, None] \
                           * u11_3[X2, :][None, None, :] * phase_m1                       
                    
                   Fc_mat_symm[N, Np, ...] = tmp1
                   
                   tmp2 = np.outer(u11_m1[X3, :], u11_m2[X2, :])[:, :, None] \
                           * u11_3[X1, :][None, None, :] * phase_1 \
                         + np.outer(u21_m1[X1, :], u11_m2[X2, :])[:, :, None] \
                           * u21_3[X3, :][None, None, :] * phase_m3 \
                         + np.outer(u11_m1[X3, :], u21_m2[X1, :])[:, :, None] \
                           * u21_3[X2, :][None, None, :] * phase_1 \
                         + np.outer(u11_m1[X2, :], u11_m2[X3, :])[:, :, None] \
                           * u11_3[X1, :][None, None, :] * phase_2 \
                         + np.outer(u11_m1[X2, :], u21_m2[X1, :])[:, :, None] \
                           * u21_3[X3, :][None, None, :] * phase_m3 \
                         + np.outer(u21_m1[X1, :], u11_m2[X3, :])[:, :, None] \
                           * u21_3[X2, :][None, None, :] * phase_2
                           
                   Fd_mat_symm[N, Np, ...] = tmp2
                   
           return Fc_mat_symm, Fd_mat_symm
        
        else:
            
           Fc_mat_symm = np.zeros((2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                               dtype=complex)
           Fd_mat_symm = np.zeros((2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                               dtype=complex)
              
           for N in range(2):
              X1 = num_sub*N + sub1
              
              for Np in range(2):
                  X2 = num_sub*Np + sub1
                  
                  for M in range(2):
                      X3 = num_sub*M + sub2
                      
                      # Fc_mat_symm is symmetric under permutation of q1, q2
                      
                      tmp1 = np.outer(u11_m1[X1, :], u21_m2[X2, :])[:, :, None] \
                           * u11_3[X3, :][None, None, :] * phase_3 \
                           + np.outer(u21_m1[X3, :], u21_m2[X2, :])[:, :, None] \
                           * u21_3[X1, :][None, None, :] * phase_m1 \
                           + np.outer(u11_m1[X1, :], u21_m2[X3, :])[:, :, None] \
                           * u11_3[X2, :][None, None, :] * phase_m2 \
                           + np.outer(u21_m1[X2, :], u11_m2[X1, :])[:, :, None] \
                           * u11_3[X3, :][None, None, :] * phase_3 \
                           + np.outer(u21_m1[X2, :], u21_m2[X3, :])[:, :, None] \
                           * u21_3[X1, :][None, None, :] * phase_m2 \
                           + np.outer(u21_m1[X3, :], u11_m2[X1, :])[:, :, None] \
                           * u11_3[X2, :][None, None, :] * phase_m1 

                         
                      Fc_mat_symm[N, Np, M, ...] = tmp1
                   
                   # Fd_mat_symm is symmetric under the permutation q1, q2                 
                   
                     
                      tmp2 = np.outer(u11_m1[X3, :], u11_m2[X2, :])[:, :, None] \
                           * u11_3[X1, :][None, None, :] * phase_1 \
                         + np.outer(u21_m1[X1, :], u11_m2[X2, :])[:, :, None] \
                           * u21_3[X3, :][None, None, :] * phase_m3 \
                         + np.outer(u11_m1[X3, :], u21_m2[X1, :])[:, :, None] \
                           * u21_3[X2, :][None, None, :] * phase_1 \
                         + np.outer(u11_m1[X2, :], u11_m2[X3, :])[:, :, None] \
                           * u11_3[X1, :][None, None, :] * phase_2 \
                         + np.outer(u11_m1[X2, :], u21_m2[X1, :])[:, :, None] \
                           * u21_3[X3, :][None, None, :] * phase_m3 \
                         + np.outer(u21_m1[X1, :], u11_m2[X3, :])[:, :, None] \
                           * u21_3[X2, :][None, None, :] * phase_2
                                                         
                      Fd_mat_symm[N, Np, M, ...] = tmp2
                      
        return Fc_mat_symm, Fd_mat_symm
                   
                   
           
    Fcmat1 = np.zeros((num_bond, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fdmat1 = np.zeros((num_bond, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fcmat2 = np.zeros((num_bond, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fdmat2 = np.zeros((num_bond, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    
    Fcmat3 = np.zeros((num_bond, 2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fdmat3 = np.zeros((num_bond, 2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fcmat4 = np.zeros((num_bond, 2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fdmat4 = np.zeros((num_bond, 2, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
        
    for bond in range(12):
        
        bond_vec = delta_ij[:, bond]
        subi = sub_idx[bond, 0]
        subj = sub_idx[bond, 1]
        
        Fcmat1[bond, :, :, :, :, :], \
            Fdmat1[bond, :, :, :, :, :] = Fcd_mat_symm(subj, subj, bond_vec, 0)
        
        Fcmat2[bond, :, :, :, :, :], \
            Fdmat2[bond, :, :, :, :, :] = Fcd_mat_symm(subi, subi, bond_vec, 0)
            
        Fcmat3[bond, :, :, :, :, :, :], \
            Fdmat3[bond, :, :, :, :, :, :] = Fcd_mat_symm(subj, subi, \
                                                          bond_vec, 2)
        
        Fcmat4[bond, :, :, :, :, :, :], \
            Fdmat4[bond, :, :, :, :, :, :] = Fcd_mat_symm(subi, subj, \
                                                          bond_vec, 1)
        
        
    tmp1 = JJJ[:, :, :, None, None, None] * Fcmat1 \
         + JJJ[:, :, :, None, None, None].conj() * Fdmat1 \
         + III[:, :, :, None, None, None] * Fcmat2 \
         + III[:, :, :, None, None, None].conj() * Fdmat2

         
    res1 = tmp1.sum(axis=(0, 1, 2))  
    
    tmp2 = JJI[:, :, :, :, None, None, None] * Fcmat3 \
         + JJI[:, :, :, :, None, None, None].conj() * Fdmat3 \
         + IIJ[:, :, :, :, None, None, None] * Fcmat4 \
         + IIJ[:, :, :, :, None, None, None].conj() * Fdmat4 

    
    res2 = tmp2.sum(axis=(0, 1, 2, 3))
    
    Fcmat5 = np.zeros((num_sub, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
    Fdmat5 = np.zeros((num_sub, 2, 2, 2*num_sub, 2*num_sub, 2*num_sub), \
                      dtype=complex)
        
    for sublat in range(num_sub):
       
        Fcmat5[sublat, :, :, :, :, :], Fdmat5[sublat, :, :, :, :, :] \
            = Fcd_mat_symm(sublat, sublat, 0.0, 0)
         
    tmp3 = III_onsite[:, :, :, None, None, None] * Fcmat5 \
        + III_onsite_c[:, :, :, None, None, None] * Fdmat5
        
    res3 = tmp3.sum(axis=(0, 1, 2))
        
    V2 = res1+2*res2+res3
    
    return V2
  

                    
        
 
                
                
                
                
        
        
        
            
        
            
        
    

        
    
