#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:23:08 2019

@author: Hao
"""

import numpy as np
import cofig as cf
import ex_info as exi
import SU3rot

num_bond = cf.num_bond
num_sub = cf.num_sub
delta_ij = cf.delta_ij
sub_idx = cf.sub_idx
f0 = cf.f0
f1 = cf.f1
f2 = cf.f2
A_mat = cf.A_mat
broadening =cf.broadening
tHij = exi.tHij
thi = exi.thi
Rmat = SU3rot.R_mat


def sw_hamiltonian(q):
    """
    the spin wave hamiltonian in momentum space
    """
    
    ham11 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    ham22 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    ham12 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    ham21 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    hos11 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    hos22 = np.zeros((2*num_sub, 2*num_sub), dtype=complex)
    
    for flavor1 in range(2):
        for flavor2 in range(2):
            
            f1m  = f1[:, flavor1]
            f1mp = f1[:, flavor2]
            f2mmp = f2[:, flavor1, flavor2]
            
            for bond in range(num_bond):
                
                bond_vec = delta_ij[:, bond]
                phase = np.exp(1j*2.0*np.pi*q @ bond_vec)
                cphase = np.conj(phase)
                sub1 = sub_idx[bond, 0]
                sub2 = sub_idx[bond, 1]
                
                THij = tHij[bond, :, :]
                
                #ham11 
                tmp = 2.0*f2mmp @ THij @ f0
                ham11[num_sub*flavor1+sub1, num_sub*flavor2+sub1] += tmp
                
                tmp = np.conj(f1m) @ THij @ f1mp * phase
                ham11[num_sub*flavor1+sub1, num_sub*flavor2+sub2] += tmp
                
                tmp = f1mp @ THij @ np.conj(f1m) * cphase
                ham11[num_sub*flavor1+sub2, num_sub*flavor2+sub1] += tmp
                
                tmp = 2.0 * f0 @ THij @ f2mmp
                ham11[num_sub*flavor1+sub2, num_sub*flavor2+sub2] += tmp
        
                #ham22
                tmp = 2.0 * np.conj(f2mmp) @ THij @ f0
                ham22[num_sub*flavor1+sub1, num_sub*flavor2+sub1] += tmp
                
                tmp = f1m @ THij @ np.conj(f1mp) * phase
                ham22[num_sub*flavor1+sub1, num_sub*flavor2+sub2] += tmp
                
                tmp = np.conj(f1mp) @ THij @ f1m * cphase
                ham22[num_sub*flavor1+sub2, num_sub*flavor2+sub1] += tmp
                
                tmp = 2.0 * f0 @ THij @ np.conj(f2mmp)
                ham22[num_sub*flavor1+sub2, num_sub*flavor2+sub2] += tmp
                
                #ham12 
                tmp = np.conj(f1m) @ THij @ np.conj(f1mp) * phase
                ham12[num_sub*flavor1+sub1, num_sub*flavor2+sub2] += tmp
                
                tmp = np.conj(f1mp) @ THij @ np.conj(f1m) * cphase
                ham12[num_sub*flavor1+sub2, num_sub*flavor2+sub1] += tmp
                
                #ham21
                tmp = f1m @ THij @ f1mp * phase
                ham21[num_sub*flavor1+sub1, num_sub*flavor2+sub2] += tmp
                
                tmp = f1mp @ THij @ f1m * cphase
                ham21[num_sub*flavor1+sub2, num_sub*flavor2+sub1] += tmp
                
            for sublat in range(num_sub):
                
                tH = thi[sublat, :]
                
                tmp = 2.0 * tH @ f2mmp
                hos11[num_sub*flavor1+sublat, num_sub*flavor2+sublat] += tmp
                
                tmp = 2.0 * tH @ np.conj(f2mmp)
                hos22[num_sub*flavor1+sublat, num_sub*flavor2+sublat] += tmp
                
    ham = np.zeros((4*num_sub, 4*num_sub), dtype=complex)          
    ham[0:2*num_sub, 0:2*num_sub] = 0.5*(ham11+hos11)
    ham[2*num_sub:4*num_sub, 2*num_sub:4*num_sub] = 0.5*(ham22+hos22)
    ham[0:2*num_sub, 2*num_sub:4*num_sub] = 0.5*ham12
    ham[2*num_sub:4*num_sub, 0:2*num_sub] = 0.5*ham21
    
    return ham


def eigensystem(q):
    
    hlsw = sw_hamiltonian(q)
    ham = A_mat @ hlsw
    eigval, eigvec = np.linalg.eig(ham)
    eigval = np.real(eigval)
    idx = eigval.argsort()[::-1]
    eigval = eigval[idx]
    eigvec = eigvec[:, idx]
    tmp = eigvec.conj().T @ A_mat @ eigvec
    
    # para-renormalization of Bogoliubov eigenvectors
    for kk in range(4*num_sub):
        eigvec[:, kk] = eigvec[:, kk]/np.sqrt(np.abs(tmp[kk, kk]))
        
    return 2.0*eigval, eigvec


def greenfunction(omega, q):
    """
    This function calculates the green function of H-P bosons
    with momentum q and energy range omega
    omega can be either a single energy or a discrete energy range,
    but it must be a numpy array
    """
    
    len_omega = len(omega)
    
    gf = np.zeros([len_omega, 4*num_sub, 4*num_sub], dtype=complex)
    ek, ubov = eigensystem(q)
    ek_m, ubov_m = eigensystem(-q)
    
    ubov_m = ubov_m.conj()
    minus_mat = np.zeros([4*num_sub, 4*num_sub], dtype=complex)
    plus_mat = np.zeros([4*num_sub, 4*num_sub], dtype=complex)
    
    u11 = ubov[:2*num_sub, :2*num_sub]
    u21 = ubov[2*num_sub:, :2*num_sub]
    u11_m = ubov_m[:2*num_sub, :2*num_sub]
    u21_m = ubov_m[2*num_sub:, :2*num_sub]
    
    for band in range(2*num_sub):
        
        minus_mat[:2*num_sub, :2*num_sub] = np.outer(u11[:, band], u11[:, band].conj())
        minus_mat[:2*num_sub, 2*num_sub:] = np.outer(u11[:, band], u21[:, band].conj())
        minus_mat[2*num_sub:, :2*num_sub] = np.outer(u21[:, band], u11[:, band].conj())
        minus_mat[2*num_sub:, 2*num_sub:] = np.outer(u21[:, band], u21[:, band].conj())
        
        plus_mat[:2*num_sub, :2*num_sub] = np.outer(u21_m[:, band], u21_m[:, band].conj())
        plus_mat[:2*num_sub, 2*num_sub:] = np.outer(u21_m[:, band], u11_m[:, band].conj())
        plus_mat[2*num_sub:, :2*num_sub] = np.outer(u11_m[:, band], u21_m[:, band].conj())
        plus_mat[2*num_sub:, 2*num_sub:] = np.outer(u11_m[:, band], u11_m[:, band].conj())
        
        tmp1 = np.reshape(omega - ek[band] + 1j*broadening, len_omega)
        tmp2 = np.reshape(omega + ek_m[band] + 1j*broadening, len_omega)
        
        # use this to avoid loop 
        gf += - minus_mat/(tmp1[:, None, None]) + plus_mat/(tmp2[:, None, None])
        
    return gf
               
# =============================================================================
#         for i in range(4*num_sub):
#             for j in range(4*num_sub):
#                 
#                 tmp1 = np.reshape(omega - ek[band] + 1j*broadening, len_omega)
#                 tmp2 = np.reshape(omega + ek_m[band] + 1j*broadening, len_omega)
#                 
#                 gf[:, i, j] += - minus_mat[i, j]/tmp1 + plus_mat[i, j]/tmp2
# =============================================================================
                
        
                

def int_mat(sub1, sub2, m, mp):
    
    """
    calculates the matrices that will be used in the intensity calculation,
    actually no need to compute this for each q, so maybe this can be optimized
    and stored in a separate file
    
    """
    
    R1t = Rmat[sub1, :, :].T.reshape(8, 8)
    R2 = Rmat[sub2, :, :].reshape(8, 8)
    r1t = R1t[:3, :]
    r2 = R2[:, :3]    
    
    f1m = f1[:, m]
    f1mp = f1[:, mp]

    mat1 = r1t @ np.outer(f1m, f1mp) @ r2
    mat2 = r1t @ np.outer(f1m.conj(), f1mp.conj()) @ r2
    mat3 = r1t @ np.outer(f1m, f1mp.conj()) @ r2
    mat4 = r1t @ np.outer(f1m.conj(), f1mp) @ r2    
    
    return mat1, mat2, mat3, mat4


def intensity(omega, qx, qy, qz):
    
    """ 
    output single-crystal neutron intensity for given qx, qy, and qz
    """
    
    qq = np.array([qx, qy, qz])
    q1, q2, q3 = cf.kxyTok12(qx, qy, qz)
    len_omega = len(omega)
    q = np.array([q1, q2, q3])
    
    chi_mat = np.zeros((len_omega, 3, 3), dtype=complex)
    
    gf = greenfunction(omega, q)
    gf11 = gf[:, :2*num_sub, :2*num_sub]
    gf12 = gf[:, :2*num_sub, 2*num_sub:]
    gf21 = gf[:, 2*num_sub:, :2*num_sub]
    gf22 = gf[:, 2*num_sub:, 2*num_sub:] 
    
    for sub1 in range(num_sub):
        for sub2 in range(num_sub):
                        
            for m in range(2):
                for mp in range(2):
                    
                    mat1, mat2, mat3, mat4 = int_mat(sub1, sub2, m, mp)
                    
                    element1 = gf12[:, 4*m+sub1, 4*mp+sub2]
                    element2 = gf21[:, 4*m+sub1, 4*mp+sub2]
                    element3 = gf11[:, 4*m+sub1, 4*mp+sub2]
                    element4 = gf22[:, 4*m+sub1, 4*mp+sub2]
                    
                    chi_mat += (-1/num_sub)\
                            * (element1[:, None, None]*mat1 + element2[:, None, None]*mat2 \
                             + element3[:, None, None]*mat3 + element4[:, None, None]*mat4) 
                             
    sqw_mat = -2.0 * np.imag(chi_mat)
    ff = cf.formfactor(qq)
    sc_inten = ff*cf.projector(qx, qy, qz, sqw_mat)
    
    return sc_inten
                     
            
    
    
    
    

            
  
            
    
    
    
    
    

    

                

               
    
                
                
                
                
                
                
                
                
            
    



