#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:05:00 2019

@author: Hao
"""

import GNLSW 
#import GLSW
import numpy as np

q0 = np.array([0.34, 0.18, 0.0])

q1 = np.array([0.45, 0.44, 0.0])
q2 = np.array([-0.11, -0.26, 0.0])
q3 = -(q1+q2)

band1 = 0
band2 = 1
band3 = 2

v1 = GNLSW.V1_cubic(band1, band2, band3, q1, q2, q3)
v2 = GNLSW.V2_cubic(band1, band2, band3, q1, q2, q3)
#eig, evec = GLSW.eigensystem(q0)
print("sink vertex=", abs(v1))
print("decay vertex=", abs(v2))

