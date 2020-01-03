#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 12:33:20 2019

@author: hao
"""

""" 
To run this program:
    python3 main_opt_Ecl_scipy.py <input.txt>,
    where <input.txt> is the input file for cofig
"""
import numpy as np
import time
import Ecl_opt_obj as Eobj
import fun_generators as fg
import cofig as cf
import scipy.optimize as sciopt


num_sub = cf.num_sub
dirpath = cf.path

x_global_opt = np.zeros(4*num_sub)
e_global_opt = 1000.0

x0 = np.zeros(4*num_sub)
lb0 = np.zeros(4*num_sub)
ub0 = np.zeros(4*num_sub)

lb0[:] = -0.001

bb = sciopt.Bounds(lb0, ub0)

for flag in range(4):
    ub0[flag*4] = 2*np.pi + 0.001
    ub0[flag*4+1] = 2*np.pi + 0.001
    ub0[flag*4+2] = np.pi + 0.001 
    ub0[flag*4+3] = 2*np.pi + 0.001
    
# number of random number experiments
N_rand = 100

for flag1 in range(N_rand):
    start_time = time.time()
    np.random.seed(seed=None)
    
    for flag in range(4):
        x0[flag*4] = 2*np.pi* np.random.rand(1)
        x0[flag*4+1] = 2*np.pi* np.random.rand(1)
        x0[flag*4+2] = np.pi* np.random.rand(1)
        x0[flag*4+3] = 2*np.pi* np.random.rand(1)
        
        
    res = sciopt.minimize(Eobj.Energy_cl, x0, method='L-BFGS-B', 
                   bounds=bb,
                   options={'ftol': 2.220446049250313e-09, 'disp': 2, 'maxfun': 20000})
    x_local_opt = res.x
    e_local_opt = res.fun
    return_flag = res.status
# =============================================================================
#     opt = nlopt.opt(nlopt.LN_NELDERMEAD, 4*num_sub)
#     opt.set_lower_bounds(lb0)
#     opt.set_upper_bounds(ub0)
#     opt.set_min_objective(Eobj.Energy_cl)
#     opt.set_xtol_rel(1e-4)
#     opt.set_maxeval = 1000
#     x_local_opt = opt.optimize(x0)
#     e_local_opt = opt.last_optimum_value()
#     return_flag = opt.last_optimize_result()
#     # nlopt return values
#     # 1: generic success; 2: stopval reached; 
#     # 3: ftol_rel (or abs) reached; 4. x_tol reached
#     # 5: maxevl reached; 6: maxtime reached
#     # negative: error, check nlopt document
# =============================================================================
    
    if e_local_opt < e_global_opt:
        print('find new local minima, updates!')
        x_global_opt = x_local_opt
        e_global_opt = e_local_opt
    
    print("finishing the %s th random experiment!" % flag1)
    print("time elpase in this experiment is %s seconds" % (time.time() - start_time))
    print("the return flag is %s" % return_flag)

mf_vals = np.zeros((4, 4))

for flag2 in range(4):
    alpha1 = x_global_opt[flag2*4]
    alpha2 = x_global_opt[flag2*4+1]
    theta  = x_global_opt[flag2*4+2]
    phi    = x_global_opt[flag2*4+3]
    
    tmp1 = fg.fun_sp(alpha1, alpha2, theta, phi)
    tmp2 = fg.fun_sm(alpha1, alpha2, theta, phi)
    mf_vals[flag2, 0] = np.real((tmp1 + tmp2)/2.0)
    mf_vals[flag2, 1] = np.real((tmp1 - tmp2)/(2.0*1j))
    tmp3 = fg.fun_sz(theta, phi)
    mf_vals[flag2, 2] = tmp3
    tmp4 = fg.fun_sz_sq(theta)
    mf_vals[flag2, 3] = tmp4
    
fname = dirpath + 'opt_angles.txt'    
np.savetxt(fname, x_global_opt)
fname1 = dirpath + 'mfvals.txt'
np.savetxt(fname1, mf_vals)  

# =============================================================================
# print("the global minimized results is: ", x_global_opt)
# print("the minimized classical energy is: ", e_global_opt)
# =============================================================================
