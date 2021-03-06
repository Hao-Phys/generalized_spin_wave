#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:05:00 2019

@author: Hao
"""

import GNLSW_vertex as vertex
import GLSW
import numpy as np
#import cofig as cf
#import ex_info as exi

path_H = 0.5
path_K = 0.0

#q0 = np.array([0.34, 0.18, 0.0])
q = np.zeros((3))
q[0] = 4.0*path_H+2.0*path_K
q[1] = path_K
q[2] = 0.5*path_K-path_H
q1 = np.array([0.23, -0.55, 0.0])
q2 = np.array([1.33, 0.87, 0.0])
q3 = -(q1+q2)

tmp, ubov1 = GLSW.eigensystem(q1)
print(tmp)
tmp, ubov2 = GLSW.eigensystem(q2)
tmp, ubov3 = GLSW.eigensystem(q3)

tmp, ubovm1 = GLSW.eigensystem(-q1)
tmp, ubovm2 = GLSW.eigensystem(-q2)
tmp, ubovm3 = GLSW.eigensystem(-q3)

num_sub = 4

band1 = 1


for band2 in range(2*num_sub):
   for band3 in range(2*num_sub):
        
       
       print("band2=", band2, "band3=", band3)
       v1 = vertex.V2_cubic_bm(band3, band2, band1, -q3, -q2, -q1, \
                            ubovm3, ubovm2, ubovm1, ubov3, ubov2, ubov1)
        #v2 = vertex.V1_cubic(band1, band2, band3, q1, q2, q3)
       #print("decay=", abs(v1))
        


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
# # =============================================================================
# #         print("real space cof is", thi[sublat, :]@ cf.f3[:, 1, mp])
# #         print("=====")
# #         print("the matrix is", U11_q3[4*mp+sublat, band3])
# # =============================================================================
#         value2 += (thi[sublat, :] @ cf.f3[:, 1, mp])*U11_q3[4*mp+sublat, band3]
#             
#     return value1, value2
# =============================================================================

#val1 = vertex.test_vertex(0, 1, 2, q1, q2, q3)
           
# =============================================================================
# print(abs(val2))
# =============================================================================
# =============================================================================
# print("source=", abs(v1))
# print("decay=", abs(v2))
# =============================================================================
