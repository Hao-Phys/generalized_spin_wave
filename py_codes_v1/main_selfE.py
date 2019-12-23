#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 10:45:08 2019

@author: hao
"""

import numpy as np
import cofig as cf
import GLSW
import GNLSW_selfE as selfE
import sys

num_sub = cf.num_sub
# input file for three momenta
inFile = sys.argv[2] 
inM = np.loadtxt(inFile)
num = inM[0]
q = np.array([inM[1], inM[2], inM[3]])

eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)

f = open('selfE.txt', 'w')
f.write('%4d' % num)
f.write('%8.3f' % q[0])
f.write('%8.3f' % q[1])
f.write('%8.3f' % q[2])

res1 = selfE.Sigma_decay(q, eq, ubov_q, ubov_mq)
res2 = selfE.Sigma_source(q, eq, ubov_q, ubov_mq)

"""
selfE.txt
ith  q1  q2  q3  e1...e8 Re(decay1) Im(decay1) ... Re(decay8) Im(decay8) ... 
Re(source1) ... Re(source8)  
"""

for flag in range(4*num_sub):
    f.write('%20.10f' %res1[flag])

for flag1 in range(2*num_sub-1):
    f.write('%20.10f' %res2[flag1])

# write change line indicator for later cat 
f.write('%20.10f\n' %res2[2*num_sub-1])
f.close()