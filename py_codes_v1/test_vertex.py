#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : test_vertex.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.14.2020
# Last Modified Date: 01.17.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:05:00 2019

@author: Hao
"""

import GNLSW_vertex as vertex
import GNLSW_vertex_new as vertex_new
import GLSW
import numpy as np
from matplotlib import pyplot as plt
# import cofig as cf
# import ex_info as exi



#q = np.array([-0.45, -0.44, 0.0])
q1 = np.array([-0.89, -0.79, 0.0])
q2 = np.array([0.11, 0.26, 0.0])
q3 = q1-q2
#q3 = np.array([0.0, 0.0, 0.0])

# q1 = np.array([0.0, 0.0, 0.0])
# q2 = np.array([0.0, 0.0, 0.0])
# q3 = q1-q2

tmp, ubov1 = GLSW.eigensystem(q1)
tmp, ubov2 = GLSW.eigensystem(q2)
tmp, ubov3 = GLSW.eigensystem(q3)

tmp, ubovm1 = GLSW.eigensystem(-q1)
tmp, ubovm2 = GLSW.eigensystem(-q2)
tmp, ubovm3 = GLSW.eigensystem(-q3)

# band1 = 1
# band2 = 5
# band3 = 7

# for band1 in range(6):
    # band2 = band1 + 2
    # band3 = 4
    # v1 = vertex.V2_cubic_bm(band1, band2, band3, q1, q2, q3, \
                        # ubov1, ubov2, ubov3, ubovm1, ubovm2, ubovm3)
# v1 = vertex.V2_cubic_bm(band3, band2, band1, -q3, -q2, q1, ubovm3, ubovm2, ubov1, \
                 # ubov3, ubov2, ubovm1)

# print(v1)

v1p = vertex_new.V_cubic_decay_bm(-q3, -q2, q1, ubovm3, ubovm2, ubov1, \
        ubov3.conj(), ubov2.conj(), ubovm1.conj())


#print(v1p[band3, band2, band1])

# print(v1p[band1, band2, band3])
res = np.zeros((64))
band1 = 4
counter = 0
for band2 in range(8):
    for band3 in range(8):
        res[counter] = abs(v1p[band3, band2, band1]) 
        # print(abs(v1p[band3, band2, band1]))
        counter += 1


data_ss = np.loadtxt('/Users/Hao/Desktop/iij_nonzeroq.txt')
kk = range(64)
plt.plot(kk, res, 'r*')
plt.plot(kk, data_ss, 'b^')
plt.show()
#print(abs(v1p[band3, band2, band1]))1)

#print(abs(v1))


# =============================================================================
# q = np.zeros((3))
# q[0] = 4.0*path_H+2.0*path_K
# q[1] = path_K
# q[2] = 0.5*path_K-path_H
# =============================================================================

# =============================================================================
# q1 = q
# q2 = np.array([1.0, 1.0, 1.0])
# q3 = -(q1+q2)
# 
# =============================================================================
# =============================================================================
# tmp, ubov1 = GLSW.eigensystem(q2)
# #print(ubov1.conj().T @ cf.A_mat @ ubov1)
# print(((ubov1[1, :]*ubov1[1, :]).conj()).sum())
# =============================================================================
#print(tmp[3])
# tmp, ubov2 = GLSW.eigensystem(q2)
# print(tmp[7])
# tmp, ubov3 = GLSW.eigensystem(q3)
# print(tmp[7])

# tmp, ubovm1 = GLSW.eigensystem(-q1)
# tmp, ubovm2 = GLSW.eigensystem(-q2)
# tmp, ubovm3 = GLSW.eigensystem(-q3)

# num_sub = 4
# # =============================================================================
# # band3 = 7
# # band2 = 7
# # =============================================================================
# band1 = 0

# =============================================================================
# v1 = vertex.V2_cubic_bm(band3, band2, band1, -q3, -q2, -q1, \
#                             ubovm3, ubovm2, ubovm1, ubov3, ubov2, ubov1)
# print(v1)
# print(abs(v1))
# =============================================================================

# for band2 in range(2*num_sub):
   # for band3 in range(2*num_sub):
        
       
       # #print("band2=", band2, "band3=", band3)
       # v1 = vertex.V2_cubic_bm(band3, band2, band1, -q3, -q2, -q1, \
                            # ubovm3, ubovm2, ubovm1, ubov3, ubov2, ubov1)
        # #v2 = vertex.V1_cubic(band1, band2, band3, q1, q2, q3)
       # if (band2 == band3):
          # print("decay=", abs(v1))
        


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
