#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:27:53 2020

@author: Hao
"""

import sys
sys.path.append('../')

import numpy as np
#import cofig as cf
import GLSW
import GNLSW_vertex_new as vertex
import GNLSW_vertex as vertex_old
from matplotlib import pyplot as plt

path_H = 0.5
path_K = 0.0

q = np.zeros((3))
q[0] = 4.0*path_H+2.0*path_K
q[1] = path_K
q[2] = 0.5*path_K-path_H

eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)

print("incoming energy is", )
kk = np.linspace(0.0, 1.0, num=40)
v1 = np.zeros((len(kk)))
v2 = np.zeros((len(kk)))

for flag in range(len(kk)):
    
    k = np.array([0.0, kk[flag], 0.0])
    q2 = k
    q3 = q - k
    e2, ubov_2 = GLSW.eigensystem(q2)
    tmp, ubov_m2 = GLSW.eigensystem(-q2)
    e3, ubov_3 = GLSW.eigensystem(q3)
    tmp, ubov_m3 = GLSW.eigensystem(-q3)
    
    V1 = vertex.V_cubic_decay(q3, q2, q, ubov_3, ubov_2, ubov_q, \
                              ubov_m3.conj(), ubov_m2.conj(), ubov_mq.conj())
    V2 = vertex.V_cubic_source(q3, q2, q, ubov_3, ubov_2, ubov_q, \
                              ubov_m3.conj(), ubov_m2.conj(), ubov_mq.conj())   
        
    v1[flag] = abs(V1[7, 7, 3])
    v2[flag] = abs(V2[7, 7, 3])
        
    
plt.plot(kk, v1, 'r-', kk, v2, 'b-')
plt.show()    
    
    

