#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:23:08 2019

@author: Hao
"""

import numpy as np
import cofig as cf
import ex_info as exi

num_bond = cf.num_bond
num_sub = cf.num_sub
delta_ij = cf.delta_ij
sub_idx = cf.sub_idx
f0 = cf.f0
f1 = cf.f1
f2 = cf.f2
tHij = exi.tHij
thi = exi.thi

# the linear spin wave Hamiltonian

def sw_hamiltonian(q):
    
    
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
                ham11[4*flavor1+sub1, 4*flavor2+sub1] += tmp
                
                tmp = np.conj(f1m) @ THij @ f1mp * phase
                ham11[4*flavor1+sub1, 4*flavor2+sub2] += tmp
                
                tmp = f1mp @ THij @ np.conj(f1m) * cphase
                ham11[4*flavor1+sub2, 4*flavor2+sub1] += tmp
                
                tmp = 2.0 * f0 @ THij @ f2mmp
                ham11[4*flavor1+sub2, 4*flavor2+sub2] += tmp
        
                #ham22
                tmp = 2.0 * np.conj(f2mmp) @ THij @ f0
                ham22[4*flavor1+sub1, 4*flavor2+sub1] += tmp
                
                tmp = f1m @ THij @ np.conj(f1mp) * phase
                ham22[4*flavor1+sub1, 4*flavor2+sub2] += tmp
                
                tmp = np.conj(f1mp) @ THij @ f1m * cphase
                ham22[4*flavor1+sub2, 4*flavor2+sub1] += tmp
                
                tmp = 2.0 * f0 @ THij @ np.conj(f2mmp)
                ham22[4*flavor1+sub2, 4*flavor2+sub2] += tmp
                
                #ham12 
                tmp = np.conj(f1m) @ THij @ np.conj(f1mp) * phase
                ham12[4*flavor1+sub1, 4*flavor2+sub2] += tmp
                
                tmp = np.conj(f1mp) @ THij @ np.conj(f1m) * cphase
                ham12[4*flavor1+sub2, 4*flavor2+sub1] += tmp
                
                #ham21
                tmp = f1m @ THij @ f1mp * phase
                ham21[4*flavor1+sub1, 4*flavor2+sub2] += tmp
                
                tmp = f1mp @ THij @ f1m * cphase
                ham21[4*flavor1+sub2, 4*flavor2+sub1] += tmp
                
            for sublat in range(num_sub):
                
                tH = thi[sublat, :]
                
                tmp = 2.0 * tH @ f2mmp
                hos11[4*flavor1+sublat, 4*flavor2+sublat] += tmp
                
                tmp = 2.0 * tH @ np.conj(f2mmp)
                hos22[4*flavor1+sublat, 4*flavor2+sublat] += tmp
                
    ham = np.zeros((4*num_sub, 4*num_sub), dtype=complex)          
    ham[0:2*num_sub, 0:2*num_sub] = 0.5*(ham11+hos11)
    ham[2*num_sub:4*num_sub, 2*num_sub:4*num_sub] = 0.5*(ham22+hos22)
    ham[0:2*num_sub, 2*num_sub:4*num_sub] = 0.5*ham12
    ham[2*num_sub:4*num_sub, 0:2*num_sub] = 0.5*ham21
    
    return ham
                
                
    
                
                
                
                
                
                
                
                
            
    



