#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : test_vertex.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.14.2020
# Last Modified Date: 02.11.2020
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
import cofig as cf
# import ex_info as exi

# qq = np.array([0.43, 0.35, 0.0])
# hswq = GLSW.sw_hamiltonian(qq)
# fname = 'tmp/hq.txt'
# np.savetxt(fname, hswq)
# hswmq = GLSW.sw_hamiltonian(-qq)
# fname = 'tmp/hmq.txt'
# np.savetxt(fname, hswmq)
# JJI = vertex.JJI
# fn1 = 'tmp/Re_real_space_vertex.txt'
# fn2 = 'tmp/Im_real_space_vertex.txt'

# JJI_bond10 = JJI[10, :, :, :]
# ReJJI_mat = np.zeros((2, 4))
# ImJJI_mat = np.zeros((2, 4))

# for boson1 in range(2):
    # ctt = 0
    # for boson2 in range(2):
        # for boson3 in range(2):
            # ReJJI_mat[boson1, ctt] = np.real(JJI_bond10[boson1, boson2, boson3])
            # ImJJI_mat[boson1, ctt] = np.imag(JJI_bond10[boson1, boson2, boson3])
            # ctt += 1

# np.savetxt(fn1, ReJJI_mat)
# np.savetxt(fn2, ImJJI_mat)
            
#q = np.array([-0.45, -0.44, 0.0])
q1 = np.array([-0.89, -0.79, 0.0])
q2 = np.array([0.11, 0.26, 0.0])
q3 = q1-q2

# q1 = np.array([0.0, 0.0, 0.0])
# q2 = np.array([0.0, 0.0, 0.0])
# q3 = q1-q2
hsw1 = GLSW.sw_hamiltonian(q1)
hsw2 = GLSW.sw_hamiltonian(q2)
hsw3 = GLSW.sw_hamiltonian(q3)

# fname = 'tmp/hamq1.txt'
# np.savetxt(fname, hsw1)

# fname = 'tmp/hamq2.txt'
# np.savetxt(fname, hsw2)

# fname = 'tmp/hamq3.txt'
# np.savetxt(fname, hsw3)

# tmp, ubov1 = GLSW.eigensystem(q1)
# fname = 'tmp/ubov_mq1.txt'
# np.savetxt(fname, ubov1)

# tmp, ubov2 = GLSW.eigensystem(q2)
# fname = 'tmp/ubov2.txt'
# np.savetxt(fname, ubov2)

# tmp, ubov3 = GLSW.eigensystem(q3)
# fname = 'tmp/ubov3.txt'
# np.savetxt(fname, ubov3)

# JJI = vertex.JJI
# print(JJI[10, 0, :, :])
# print(JJI[10, 1, :, :])

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

# v1p = vertex_new.V_cubic_decay_bm(-q3, -q2, q1, ubovm3, ubovm2, ubov1, \
        # ubov3.conj(), ubov2.conj(), ubovm1.conj())

# v1p = vertex_new.V_cubic_decay_bm(q3, q2, -q1, ubov3, ubov2, ubovm1, \
        # ubovm3.conj(), ubovm2.conj(), ubov1.conj())

# print(v1p[band3, band2, band1])

# band1 = 4
# band2 = 0
# band3 = 1

# vec = np.array([0.25, 0.0, 0.0])
# va1, va2, va3, va = vertex.Fd_fun1(0, 0, 3, band1, band2, band3, \
                      # ubov1, ubov2, ubov3, -q1, q2, q3, vec, 2)

# print('va1=', va1)
# print('va2=', va2)
# print('va3=', va3)
# print('va=', va)

res = np.zeros((64))
res1 = np.zeros((64))

v1 = vertex_new.V_cubic_decay_bm(-q3, -q2, q1, ubovm3, ubovm2, ubov1, \
        ubov3.conj(), ubov2.conj(), ubovm1.conj())

band1 = 4
counter = 0
for band2 in range(8):
    for band3 in range(8):

        # v1 = vertex.V2_cubic_bm(band3, band2, band1, -q3, -q2, q1, \
                                # ubovm3, ubovm2, ubov1, ubov3, ubov2, ubovm1)
        res[counter] = abs(v1[band3, band2, band1])
        # res1[counter] = np.imag(v1)
        counter += 1


#data_ss = np.loadtxt('/Users/Hao/Desktop/c3b.txt')
# data_ss = np.loadtxt('/Users/Hao/Desktop/vertex_jji_sublattice1_bond2_band1=5.txt')
fname1 = '/Users/Hao/Desktop/all_on1.txt'
np.savetxt(fname1, res)
# fname2 = '/Users/Hao/Desktop/vvv/Imbond2.txt'
# np.savetxt(fname2, res1)
# kk = range(64)
# plt.plot(kk, res, 'r*')
# plt.plot(kk, data_ss[:, 10], 'b^')
# plt.show()
# print(abs(v1p[band3, band2, band1]))1)

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
