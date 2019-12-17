#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 20:36:07 2019

@author: hao
"""

import numpy as np
import GLSW
import GNLSW_vertex as vertex
import cofig as cf
import time

num_sub = cf.num_sub
q = np.array([0.12, 0.31, 0.0])
eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)
omega = 0.12

k = np.array([0.23,0.54,0.0])
qmk = q - k

st = time.time()
result = 0.0
ek, ubov_k = GLSW.eigensystem(k)
eqmk, ubov_qmk = GLSW.eigensystem(qmk)
tmp, ubov_mk = GLSW.eigensystem(-k)
tmp, ubov_mqmk = GLSW.eigensystem(-qmk)


for band2 in range(2*num_sub):
    for band3 in range(2*num_sub):
        
        vd = vertex.V2_cubic(0, band2, band3, q, k, qmk, \
                             ubov_q, ubov_k, ubov_qmk, \
                             ubov_mq, ubov_mk, ubov_mqmk)
            
        tmp = vd.conj() * vd/(omega - ek[band2] - eqmk[band3] + 1j*cf.convergence)
        result += tmp
        
et = time.time()        
print(result)
print('evaluate integrand once takes time = ', et-st, 's')