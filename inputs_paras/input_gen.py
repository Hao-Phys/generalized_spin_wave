#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : input_gen.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.13.2020
# Last Modified Date: 02.13.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 13:39:46 2019

@author: hao
"""

import os
import numpy as np

field = np.array([0.0, 1.0, 3.0, 4.0, 4.74])
#os.chdir('gs_info/')
opt_file = np.loadtxt('fit_results.txt')

J1pm   = opt_file[0]
J1zz   = opt_file[1]
J1pmpm = opt_file[2]
J1zpm  = opt_file[3]
J1 = J1pm
Delta1 = J1zz/J1

J2pm   = opt_file[4]
J2zz   = opt_file[5]
J2 = J2pm
Delta2 = J2zz/J2

J3pm   = opt_file[6]
J3zz   = opt_file[7]
J3pmpm = opt_file[8]
J3zpm  = opt_file[9]
J3 = J3pm
Delta3 = J3zz/J3

Jp0pm  = opt_file[10]
Jp0zz  = opt_file[11]
Deltap0 = Jp0zz/Jp0pm

Jp1pm = opt_file[12]
Jp1zz = opt_file[13]
Deltap1 = Jp1zz/Jp1pm

Jp2apm = opt_file[14]
Jp2azz = opt_file[15]
Deltap2a = Jp2azz/Jp2apm

Jp2bpm = opt_file[16]
Jp2bzz = opt_file[17]
Deltap2b = 0.0

D_ion = opt_file[18]


for h in field:
    dir = 'h_' + str(h) + 'T'
    if os.path.isdir(dir):
        pass
    else:
        os.mkdir(dir)
        
    os.chdir(dir)
    fname = 'input1.txt'
    # if os.path.exists(fname):
        # pass
    # else:
        # paras = np.array([-2.77, 0.8884, 0.272, 2.6434, 1.955, 2.4156, 
                          # 25.729, -2.02, 0.0, -3.017, 0.0, 0.3340, 0.0,
                          # 0.1700, 0.0, 0.72, 0.1986, 0.0, 0.0, h])
    paras = np.array([J1, Delta1, J2, Delta2, J3,Delta3, 
                      D_ion, J1pm, J3pm, J1zpm, J3zpm, Jp0pm, Deltap0,
                      Jp1pm, Deltap1, Jp2apm, Deltap2a, Jp2bpm, Deltap2b, h])

    np.savetxt(fname, paras)
    
    os.chdir('../')
    
os.chdir('../')



    
    
