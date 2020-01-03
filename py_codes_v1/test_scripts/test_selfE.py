#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:16:55 2019

@author: Hao
"""

import numpy as np
import GLSW
import GNLSW_selfE as selfE
import time

st = time.time()
q = np.array([0.1, 0.2, 0.0])
eq, ubov_q = GLSW.eigensystem(q)
tmp, ubov_mq = GLSW.eigensystem(-q)

E1 = selfE.Sigma_decay(q, eq, ubov_q, ubov_mq)
et = time.time()
print(E1)
print('time elapse =', et-st, 's')