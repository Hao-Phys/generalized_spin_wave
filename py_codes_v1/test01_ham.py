#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:43:53 2019

@author: Hao
"""

import GLSW
import numpy as np
import cofig as co

qq = np.array([0.0, 0.0, 0.0])
A_mat = co.A_mat
hamm = GLSW.sw_hamiltonian(qq)
eigval, eigvec = np.linalg.eig(A_mat@hamm)
idx = eigval.argsort()[::-1]
eigval = eigval[idx]
eigvec = eigvec[idx]
print(eigval)