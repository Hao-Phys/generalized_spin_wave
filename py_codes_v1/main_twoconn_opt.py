#!/usr/bin/env python3
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
fileName = 'twoconn_h' + str(cf.field) + 'T.txt'
f = open(fileName, 'w')
f.write('%4d' % num)
f.write('%8.3f' % qq[0])
f.write('%8.3f' % qq[1])
f.write('%8.3f' % qq[2])

def two_conn(k, q, band1, band2):
    
    qmk = q - k
    ek, tmp = GLSW.eigensystem(k)
    eqmk, tmp = GLSW.eigensystem(qmk)
    funval = ek[band1] + eqmk[band2]
    
    return funval

k_global_opt = np.zeros(3)
tc_global_opt = 1000.0

k0 = np.zeros(3)
lb0 = np.array([-0.1, -0.1, -0.1])
ub0 = np.array([1.1, 1.1, 1.1])
bb = sciopt.Bounds(lb0, ub0)

N_rand = 1

for band1 in range(2*num_sub):
    for band2 in range(2*num_sub):

        print('band1=', band1)
        print('band2=', band2)
        ist = time.time()

        for flag1 in range(N_rand):
            
            np.random.seed(seed=None)
            
            k0[0] = np.random.rand(1)
            k0[1] = np.random.rand(1)
            k0[2] = np.random.rand(1)
            
            res = sciopt.minimize(two_conn, k0, args=(qq, band1, band2), 
                                  method='L-BFGS-B', 
                                  bounds=bb,
                                  options={'ftol': 2.220446049250313e-09, 
                                           'disp': None, 'maxfun': 20000})
            k_local_opt = res.x
            tc_local_opt = res.fun
            return_flag = res.status
            print('the return flag is %s' %return_flag)
            if tc_local_opt < tc_global_opt:
                print('find new local minima, update')
                k_global_opt = k_local_opt
                tc_global_opt = tc_local_opt
                
        if (band1 == 2*num_sub-1 and  band2 == 2*num_sub-1):
            f.write('%20.10f\n' %tc_global_opt)
        else:
            f.write('%20.10f' %tc_global_opt)
        iet = time.time()
        print('find one minimal, time is', iet-ist)
            

print('time elapse is %s seconds' %(time.time() - start_time))  
  
        

    

    

    
    
    

