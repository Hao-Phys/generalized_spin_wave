#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:43:53 2019

@author: Hao
"""

import GLSW
import numpy as np
import cofig as cf

qq = np.array([0.2, 0.12, 0.0])
hsw = GLSW.sw_hamiltonian(qq)
omega = np.array([0.13])
gf1 = -(omega + 1j*cf.broadening) * cf.A_mat + 2.0*hsw
gf2 = GLSW.greenfunction(omega, qq)
#eigval, eigvec = GLSW.eigensystem(qq)
resd = gf1 @ gf2[0, :, :]
print(resd - np.eye(16))