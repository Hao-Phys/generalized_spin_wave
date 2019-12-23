#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 16:56:40 2019

@author: hao
"""
import numpy as np
import cofig as cf


def R_O2O_mat(alpha1, alpha2, theta, phi):
    """ 
    8*8 matrix transforms the 8 generators from 
    global to local reference frame
    """
    
    mat = np.zeros((8, 8))
    
    fp = 1/np.sqrt(2)*(np.cos(theta) + np.sin(theta))
    fm = 1/np.sqrt(2)*(np.cos(theta) - np.sin(theta))

    F1p = 1/np.sqrt(2)*(np.cos(phi)*np.cos(alpha1) \
                    + np.sin(phi)*np.cos(alpha2))
    F1m = 1/np.sqrt(2)*(np.cos(phi)*np.cos(alpha1) \
                    - np.sin(phi)*np.cos(alpha2))

    F2p = 1/np.sqrt(2)*(np.cos(phi)*np.sin(alpha1) \
                    + np.sin(phi)*np.sin(alpha2))
    F2m = 1/np.sqrt(2)*(np.cos(phi)*np.sin(alpha1) \
                    - np.sin(phi)*np.sin(alpha2))

    F3p = (np.cos(phi))**2*np.cos(2*alpha1) \
          + (np.sin(phi))**2*np.cos(2*alpha2)
    F3m = (np.cos(phi))**2*np.cos(2*alpha1) \
          - (np.sin(phi))**2*np.cos(2*alpha2)

    F4p = (np.cos(phi))**2*np.sin(2*alpha1) \
          + (np.sin(phi))**2*np.sin(2*alpha2)
    F4m = (np.cos(phi))**2*np.sin(2*alpha1) \
          - (np.sin(phi))**2*np.sin(2*alpha2)

    mat[0, :] = np.array([fm*F1m, -fm*F2p, -fp*np.cos(alpha1+alpha2)*np.sin(2*phi), \
                        -fm*F2m, fm*F1p, fp*F3m, -fp*F4p, 0])
             
    mat[1, :] = np.array([fp*F2m, fp*F1p, fm*np.sin(alpha1+alpha2)*np.sin(2*phi), \
                        fp*F1m, fp*F2p, -fm*F4m, -fm*F3p, 0])
    
    mat[2, :] = np.array([np.sin(2*theta)*F1p, -np.sin(2*theta)*F2m, \
                        -0.5*np.cos(2*theta)*np.cos(2*phi), np.sin(2*theta)*F2p, \
                        -np.sin(2*theta)*F1m, \
                        -0.5*np.cos(2*theta)*np.cos(alpha1-alpha2)*np.sin(2*phi), \
                         0.5*np.cos(2*theta)*np.sin(alpha1-alpha2)*np.sin(2*phi), \
                         -0.5*np.sqrt(3)*np.cos(2*theta)])
    
    mat[3, :] = np.array([-fm*F2m, -fm*F1p, fp*np.sin(alpha1+alpha2)*np.sin(2*phi), \
                        -fm*F1m, -fm*F2p, -fp*F4m, -fp*F3p, 0])
    
    mat[4, :] = np.array([-fp*F1m, fp*F2p, -fm*np.cos(alpha1+alpha2)*np.sin(2*phi), \
                        fp*F2m, -fp*F1p, fm*F3m, -fm*F4p, 0])
    
    mat[5, :] = np.array([np.cos(2*theta)*F1p, -np.cos(2*theta)*F2m, \
                         0.5*np.sin(2*theta)*np.cos(2*phi), np.cos(2*theta)*F2p, \
                         -np.cos(2*theta)*F1m, \
                         0.5*np.sin(2*theta)*np.cos(alpha1-alpha2)*np.sin(2*phi), \
                        -0.5*np.sin(2*theta)*np.sin(alpha1-alpha2)*np.sin(2*phi), \
                         0.5*np.sqrt(3)*np.sin(2*theta)])
    
    mat[6, :] = np.array([-F2p, -F1m, 0, F1p, F2m, 0, 0, 0])
    
    mat[7, :] = np.array([0, 0, 0.5*np.sqrt(3)*np.cos(2*phi), 0, 0, \
                         0.5*np.sqrt(3)*np.cos(alpha1-alpha2)*np.sin(2*phi), \
                        -0.5*np.sqrt(3)*np.sin(alpha1-alpha2)*np.sin(2*phi), -0.5])
    return mat

""" 
SU(3) global to local frame rotation
"""


fname = cf.path + 'opt_angles.txt'
angles = np.loadtxt(fname)
num_sub = cf.num_sub

R_mat = np.zeros((num_sub, 8, 8))

for ii in range(num_sub):
    alpha1 = angles[4*ii]
    alpha2 = angles[4*ii+1]
    theta  = angles[4*ii+2]
    phi    = angles[4*ii+3]
    R_mat[ii, :, :] = R_O2O_mat(alpha1, alpha2, theta, phi)