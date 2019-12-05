#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 13:39:46 2019

@author: hao
"""

import os
import numpy as np

field = np.array([0.0, 1.0, 3.0, 4.0, 4.74])
os.chdir('gs_info/')

for h in field:
    dir = 'h=' + str(h) + 'T'
    if os.path.isdir(dir):
        pass
    else:
        os.mkdir(dir)
        
    os.chdir(dir)
    fname = 'input.txt'
    if os.path.exists(fname):
        pass
    else:
        paras = np.array([-2.77, 0.8884, 0.272, 2.6434, 1.955, 2.4156, 
                          25.729, -2.02, 0.0, -3.017, 0.0, 0.3340, 0.0,
                          0.1700, 0.0, 0.72, 0.1986, 0.0, 0.0, h])
        np.savetxt(fname, paras)
    
    os.chdir('../')
    
os.chdir('../')



    
    
