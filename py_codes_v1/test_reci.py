#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : test_reci.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.17.2020
# Last Modified Date: 02.17.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import GLSW

num_sub = 4
q = np.array([0.1, 0.1, 0.1]) 
G1 = np.array([1.0, 0.0, 0.0])
G2 = np.array([2.0, 1.0, 3.0])
G3 = np.array([0.0, 0.0, 1.0])

q1 = q + G1 
q2 = q + G2
q3 = q + G3

e1, ubov1 = GLSW.eigensystem(q1)
e2, ubov2 = GLSW.eigensystem(q2)
e3, ubov3 = GLSW.eigensystem(q3)

ubov1_11 = ubov1[:2*num_sub, :2*num_sub]
ubov2_11 = ubov2[:2*num_sub, :2*num_sub]
ubov1_21 = ubov1[2*num_sub:, :2*num_sub]
ubov2_21 = ubov2[2*num_sub:, :2*num_sub]

print('ubov11/ubov11')
print(ubov1_11[:, 0]/ubov2_11[:, 0])
print('difference between e2, e3 is ')
print(e2 - e3)
