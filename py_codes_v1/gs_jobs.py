#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : gs_jobs.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.13.2020
# Last Modified Date: 02.13.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 15:36:05 2019

@author: hao
"""

import os
import numpy as np
import multiprocessing as mp

h = np.array([0.0, 1.0, 3.0, 4.0, 4.74])
lenh = len(h)
xs = range(lenh)

def sub_jobs(ff):
    hh = h[ff]
    inFile = '../inputs_paras/h_' + str(hh) + 'T/input1.txt'
    outFile = '../inputs_paras/h_' + str(hh) + 'T/output1.txt'
    cmd = 'python3 main_opt_Ecl_scipy.py ' + inFile +'>' + outFile
    os.system(cmd)
 
# =============================================================================
# hh = h[0]
# inFile = 'gs_info/h=' + str(hh) + 'T/input.npy'
# #outFile = 'gs_info/h=' + str(hh) + 'T/output.txt'
# cmd = 'python3 main_opt_Ecl.py ' + inFile 
# os.system(cmd)
# =============================================================================
  
pool = mp.Pool(processes=3)

pool.map(sub_jobs, xs)
    
pool.close()    

