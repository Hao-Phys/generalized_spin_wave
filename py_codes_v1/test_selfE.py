#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:16:55 2019

@author: Hao
"""

import numpy as np
import GNLSW_selfE as selfE

q1 = np.array([0.0, 0.0, 0.0])
band = 0
omega = 0

E1 = selfE.Sigma_decay(q1, omega, band)
print(E1)