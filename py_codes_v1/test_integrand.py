#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 20:36:07 2019

@author: hao
"""

import numpy as np
import GLSW
import GNLSW_vertex_new as vertex
import GNLSW_vertex as vertex_old
import cofig as cf
import time


num_sub = cf.num_sub
st = time.time()

omega = 0.12
q = np.array([0.12, 0.53, 0.0])
eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)

k = np.array([0.2, 0.1, 0.0])
qmk = q - k

ek, ubov_k = GLSW.eigensystem(k)
eqmk, ubov_qmk = GLSW.eigensystem(qmk)
tmp, ubov_mk = GLSW.eigensystem(-k)
tmp, ubov_mqmk = GLSW.eigensystem(-qmk)

V2_old = vertex_old.V2_cubic_bm(5, 6, 6, q, k, qmk, \
                             ubov_q, ubov_k, ubov_qmk, \
                             ubov_mq, ubov_mk, ubov_mqmk)

V2 = vertex.V_cubic_decay(q, k, qmk, ubov_q, ubov_k, ubov_qmk, \
                     ubov_mq.conj(), ubov_mk.conj(), ubov_mqmk.conj())
    
    
print('the old vertex =', V2_old)
print('the new vertex =', V2[5, 6, 6])
result = 0.0


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
        
et = time.time()        
# print(result)
print('evaluate integrand once takes time = ', et-st, 's')