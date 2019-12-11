#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:05:00 2019

@author: Hao
"""

import GNLSW_vertex
#import GLSW
import numpy as np
#import cofig as cf
import ex_info as exi


q0 = np.array([0.34, 0.18, 0.0])
q1 = np.array([0.0, 0.0, 0.0])
q2 = np.array([0.0, 0.0, 0.0])
q3 = -(q1+q2)
num_sub = 4
band1 = 0
band2 = 1
band3 = 2
thi = exi.thi


v1 = GNLSW_vertex.V1_cubic(band1, band2, band3, q1, q2, q3)
v2 = GNLSW_vertex.V2_cubic(band1, band2, band3, q1, q2, q3)
#eig, evec = GLSW.eigensystem(q0)



# =============================================================================
# def testvert(q1, q2, q3, band1, band2, band3, sublat):
#     
#     value1 = 0
#     value2 = 0
# #    ham1 = GLSW.sw_hamiltonian(-q0)
#     #print(ham-ham1)
#     
#     #ham = GLSW.sw_hamiltonian(q1)
#     eig1, vec1 = GLSW.eigensystem(q1)
#     eig2, vec2 = GLSW.eigensystem(-q1)
#     
#     eig3, vec3 = GLSW.eigensystem(q2)
# 
#     eig4, vec4 = GLSW.eigensystem(-q2)
#                               
#     eig5, vec5 = GLSW.eigensystem(q3)
#     eig6, vec6 = GLSW.eigensystem(-q3)
#     
#     #print(vec1.conj().T@cf.A_mat@vec1)
#     U11_mq1 = vec2[:2*num_sub, :2*num_sub]
#     #print(U11_mq1[:, 0].conj())
#     U11_q2 = vec3[:2*num_sub, :2*num_sub]
#     #print("====")
#     #print(U11_q2[:, 0])
#     U11_q3 = vec5[:2*num_sub, :2*num_sub]
#     
#     for m in range(2):
# 
#         value1 += U11_mq1[4*m+sublat, band1].conj()*U11_q2[4*m+sublat, band2]
#         
#     for mp in range(2):
#         print("real space cof is", thi[sublat, :]@ cf.f3[:, 1, mp])
#         print("=====")
#         print("the matrix is", U11_q3[4*mp+sublat, band3])
#         value2 += (thi[sublat, :] @ cf.f3[:, 1, mp])*U11_q3[4*mp+sublat, band3]
#             
#     return value1, value2
# =============================================================================

#val1, val2 = testvert(q1, q2, q3, 0, 1, 2, 3)
print("source=", abs(v1))
print("decay=", abs(v2))
