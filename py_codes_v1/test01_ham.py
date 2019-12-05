#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:43:53 2019

@author: Hao
"""

import GLSW
import numpy as np
#import cofig as co

qq = np.array([0.2, 0.12, 0.0])
eigval, eigvec = GLSW.eigensystem(qq)
print(eigval)