#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:05:00 2019

@author: Hao
"""

import GNLSW 
import GLSW
import numpy as np

q1 = np.array([0.0, 0.0, 0.0])
q2 = np.array([0.0, 0.0, 0.0])
q3 = -(q1+q2)

band1 = 0
band2 = 1
band3 = 2
value = GNLSW.V1_cubic(band1, band2, band3, q1, q2, q3)
eig, evec = GLSW.eigensystem(q1)

