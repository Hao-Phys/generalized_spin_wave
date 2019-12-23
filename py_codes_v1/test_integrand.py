#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 20:36:07 2019

@author: hao
"""

import numpy as np
#import numpy.matlib
#import numpy.testing as npt
import GLSW
import GNLSW_vertex_new as vertex
import GNLSW_vertex as vertex_old
import cofig as cf
import time



num_sub = cf.num_sub

# the self energy input
omega = np.array([0.12, 0.22, 0.13])
# incoming particle
q = np.array([0.12, 0.53, 0.0])
eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)
# intermediate particles
k = np.array([0.2, 0.1, 0.0])
qmk = q - k
ek, ubov_k = GLSW.eigensystem(k)
eqmk, ubov_qmk = GLSW.eigensystem(qmk)
tmp, ubov_mk = GLSW.eigensystem(-k)
tmp, ubov_mqmk = GLSW.eigensystem(-qmk)

# =============================================================================
# V2_old = vertex_old.V2_cubic_bm(1, 7, 6, q, k, qmk, \
#                              ubov_q, ubov_k, ubov_qmk, \
#                              ubov_mq, ubov_mk, ubov_mqmk)
# =============================================================================

"""
numpy-like integrand
"""
st1 = time.time()
V2 = vertex.V_cubic_source(q, k, qmk, ubov_q, ubov_k, ubov_qmk, \
                     ubov_mq.conj(), ubov_mk.conj(), ubov_mqmk.conj())
    
V2_old = vertex_old.V1_cubic(2, 3, 7, q, k, qmk, ubov_q, ubov_k, ubov_qmk, \
                     ubov_mq, ubov_mk, ubov_mqmk)   
    
print('new vertex = ', V2[2, 3, 7])
print('old vertex = ', V2_old)

# the denominator
# =============================================================================
# tmpmat1 = ek[:2*num_sub][:, None]
# tmpmat2 = eqmk[:2*num_sub][None, :]
# #tmpmat1 = np.matlib.repmat(ek[:2*num_sub], 2*num_sub, 1).T
# #tmpmat2 = np.matlib.repmat(eqmk[:2*num_sub], 2*num_sub, 1)
# esum = tmpmat1 + tmpmat2
# denomin = omega[:, None, None, None] - esum[None, None, :, :]
# numer = V2[None, :, :, :]
# tmp4D = numer.conj() * numer /denomin
# IntG = tmp4D.sum(axis=(2, 3))
# et1 = time.time()
# print('takes time', et1-st1)
# =============================================================================

"""
loop integrand
"""
# =============================================================================
# st2 = time.time()
# IntG_loop = np.zeros((len(omega), 2*num_sub), dtype=complex)
# 
# for ome in range(len(omega)):
#     for band1 in range(2*num_sub):
#         result = 0.0
#         
#         for band2 in range(2*num_sub):
#             for band3 in range(2*num_sub):
#                 vd = vertex_old.V2_cubic_bm(band1, band2, band3, q, k, qmk, \
#                                             ubov_q, ubov_k, ubov_qmk, \
#                                             ubov_mq, ubov_mk, ubov_mqmk)
#                 tmp = vd.conj() * vd/(omega[ome] - ek[band2] - eqmk[band3])
#                 result += tmp
#                 
#         IntG_loop[ome, band1] = result
#         
# et2 = time.time()
#                
#                 
# npt.assert_almost_equal(IntG, IntG_loop)
# 
# 
# print('the numpy way is', (et2-st2)/(et1-st1), 'faster than the loop way')
# =============================================================================

    
# =============================================================================
# print('the old vertex =', V2_old)
# print('the new vertex =', V2[1, 7, 6])
# =============================================================================




# =============================================================================
# for band2 in range(2*num_sub):
#     for band3 in range(2*num_sub):
#         
#         vd = vertex_old.V2_cubic_bm(0, band2, band3, q, k, qmk, \
#                              ubov_q, ubov_k, ubov_qmk, \
#                              ubov_mq, ubov_mk, ubov_mqmk)
#             
#         tmp = vd.conj() * vd/(omega - ek[band2] - eqmk[band3] + 1j*cf.convergence)
#         result += tmp
# =============================================================================
        
#et = time.time()        
# print(result)
#print('evaluate integrand once takes time = ', et-st, 's')