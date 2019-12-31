#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 22:10:40 2019

@author: hao
"""

import numpy as np 
import os

steps = 100

counter = 0.0
K = np.linspace(-1.0, 1.0, 21)
H = np.linspace(-2.0, 2.0, steps)

for flag1 in range(len(K)):
    dirK = 'K_' + str('%2.1f'%K[flag1])
    if (os.path.exists(dirK)):
        os.chdir(dirK)
    else:
        os.mkdir(dirK)
        os.chdir(dirK)
    
    counter = 0

    for flag2 in range(len(H)):
        counter += 1
        dirJob = 'job' + str(counter)
        if (os.path.exists(dirJob)):
            os.chdir(dirJob)
        else:
            os.mkdir(dirJob)
            os.chdir(dirJob)
        
        
        k1 = 2.0*H[flag2] + 4.0*K[flag1]
        k2 = H[flag2]
        k3 = 0.5*H[flag2] - K[flag1]
        
        fname = open('inMo1.txt', 'w')
        fname.write('%4d\n' %counter)
        fname.write('%20.10f\n' %k1)
        fname.write('%20.10f\n' %k2)
        fname.write('%20.10f\n' %k3)
        fname.close()
        os.chdir('../')
        
    os.chdir('../')
    
    
        
