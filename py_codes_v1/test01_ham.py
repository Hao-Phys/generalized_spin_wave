#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:43:53 2019

@author: Hao
"""

import GLSW
import numpy as np
import cofig as cf

qq = np.array([0.1, 0.2, 0.3])
hsw = GLSW.sw_hamiltonian(qq)
omega = np.array([4.427])
#gf1 = -(omega + 1j*cf.broadening) * cf.A_mat + 2.0*hsw
qx, qy, qz = cf.k12Tokxy(qq[0], qq[1], qq[2]) 
chi, sqw = GLSW.intensity(omega, qx, qy, qz)
#print(chi)
print(sqw)
#eigval, eigvec = GLSW.eigensystem(qq)
#resd = gf1 @ gf2[0, :, :]
#print(resd - np.eye(16))