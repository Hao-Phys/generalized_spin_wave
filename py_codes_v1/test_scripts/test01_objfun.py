#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 23:37:05 2019

@author: hao
"""

import Ecl_opt_obj as obj
import numpy as np
 

#infile = sys.argv[1]
vec = np.loadtxt('opt_input.txt')
vec1 = np.loadtxt('gs_info/h=0.0T/opt_angles.txt')
funval = obj.Energy_cl(vec)
funval1 = obj.Energy_cl(vec1)
print("the optimized energy from matlab= ", funval)
print("the optimized energy from scipy= ", funval1)
 