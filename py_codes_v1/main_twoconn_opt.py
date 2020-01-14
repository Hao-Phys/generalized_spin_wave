#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : main_twoconn_opt.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.11.2020
# Last Modified Date: 01.11.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:11:28 2020

@author: Hao
"""

import numpy as np
import GLSW
import time
import scipy.optimize as sciopt
import sys
import cofig as cf


inFile = sys.argv[2]
inM = np.loadtxt(inFile)
num = inM[0]
qq = np.array([inM[1], inM[2], inM[3]])

start_time = time.time()
num_sub = cf.num_sub
fileName1 = 'min_twoconn_h' + str(cf.field) + 'T.txt'
fileName2 = 'max_twoconn_h' + str(cf.field) + 'T.txt'

f1 = open(fileName1, 'w')
f1.write('%4d' % num)
f1.write('%8.3f' % qq[0])
f1.write('%8.3f' % qq[1])
f1.write('%8.3f' % qq[2])

f2 = open(fileName2, 'w')
f2.write('%4d' % num)
f2.write('%8.3f' % qq[0])
f2.write('%8.3f' % qq[1])
f2.write('%8.3f' % qq[2])

def two_conn_min(k, q, band1, band2):
    
    qmk = q - k
    ek, tmp = GLSW.eigensystem(k)
    eqmk, tmp = GLSW.eigensystem(qmk)
    funval = ek[band1] + eqmk[band2]
    
    return funval

def two_conn_max(k, q, band1, band2):
    
    qmk = q - k
    ek, tmp = GLSW.eigensystem(k)
    eqmk, tmp = GLSW.eigensystem(qmk)
    funval = -(ek[band1] + eqmk[band2])
    
    return funval

k_global_min = np.zeros(3)
k_global_max = np.zeros(3)


k0 = np.zeros(3)
lb0 = np.array([-0.1, -0.1, -0.1])
ub0 = np.array([1.1, 1.1, 1.1])
bb = sciopt.Bounds(lb0, ub0)

N_rand = 15

for band1 in range(2*num_sub):
    for band2 in range(2*num_sub):

        print('band1=', band1)
        print('band2=', band2)
        ist = time.time()
        tc_global_min = 1000.0

        for flag1 in range(N_rand):
            
            np.random.seed(seed=None)
            
            k0[0] = np.random.rand(1)
            k0[1] = np.random.rand(1)
            k0[2] = np.random.rand(1)
            
            res = sciopt.minimize(two_conn_min, k0, args=(qq, band1, band2), 
                                  method='L-BFGS-B', 
                                  bounds=bb,
                                  options={'ftol': 2.220446049250313e-09, 
                                           'disp': None, 'maxfun': 20000})
            k_local_min = res.x
            tc_local_min = res.fun
            return_flag = res.status
            print('the return flag is %s' %return_flag)
            if tc_local_min < tc_global_min:
                print('find new local minima, update')
                k_global_min = k_local_min
                tc_global_min = tc_local_min
                
        if (band1 == 2*num_sub-1 and  band2 == 2*num_sub-1):
            f1.write('%20.10f\n' %tc_global_min)
        else:
            f1.write('%20.10f' %tc_global_min)
        iet = time.time()
        print('find one minimal, time is', iet-ist)
            

for band1 in range(2*num_sub):
    for band2 in range(2*num_sub):

        print('band1=', band1)
        print('band2=', band2)
        ist = time.time()
        tc_global_max = 1000.0 

        for flag1 in range(N_rand):
            
            np.random.seed(seed=None)
            
            k0[0] = np.random.rand(1)
            k0[1] = np.random.rand(1)
            k0[2] = np.random.rand(1)
            
            res = sciopt.minimize(two_conn_max, k0, args=(qq, band1, band2), 
                                  method='L-BFGS-B', 
                                  bounds=bb,
                                  options={'ftol': 2.220446049250313e-09, 
                                           'disp': None, 'maxfun': 20000})
            k_local_max = res.x
            tc_local_max = res.fun
            return_flag = res.status
            print('the return flag is %s' %return_flag)
            if tc_local_max < tc_global_max:
                print('find new local minima, update')
                k_global_max = k_local_max
                tc_global_max = tc_local_max
                
        if (band1 == 2*num_sub-1 and  band2 == 2*num_sub-1):
            f2.write('%20.10f\n' %(-tc_global_max))
        else:
            f2.write('%20.10f' %(-tc_global_max))
        iet = time.time()
        print('find one maximal, time is', iet-ist)

print('time elapse is %s seconds' %(time.time() - start_time))  

f1.close()
f2.close()
  
        

    

    

    
    
    

