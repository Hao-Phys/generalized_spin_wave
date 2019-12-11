#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:27:15 2019

@author: Hao
"""

import numpy as np
import cofig as cf
import GLSW
import GNLSW_vertex as vertex
import pycuba

cgf = cf.convergence
num_sub = cf.num_sub


# print function copied from cuba_domo

def print_header(name):
  print('-------------------- %s test -------------------' % name)
def print_results(name, results):
  keys = ['nregions', 'neval', 'fail']
  text = ["%s %d" % (k, results[k]) for k in keys if k in results]
  print("%s RESULT:\t" % name.upper() + "\t".join(text))
  for comp in results['results']:
    print("%s RESULT:\t" % name.upper() + \
	"%(integral).8f +- %(error).8f\tp = %(prob).3f\n" % comp)
    
### specific for cuba ###

"""
vegas:
    pycuba.Vegas(integrand, ndim, userdata, epsrel=0.001, epsabs=1e-12, \
                 verbose=0, ncomp=1, seed=None, minevel=0, maxeval=50000, \
                 nstart=1000, nincrease=500, nbatch=1000, gridno=0, statefile=, \
                 nvec=1)
"""
# cuba variables

NDIM = 3
EPSREAL = 1e-3
EPSABS = 1e-12
NCOMP = 1
NCOMP2 = 2
MINEVEL = 0
MAXEVAL = 50000


def Sigma_decay(q, omega, band1):
    """
    calculates the self energy from the decay channel
    """
    
    # define the integrand as a nested function
    def integrand_fun1(ndim, xx, ncopm, ff, userdata):
        """
        the decay integrand
        """
        
        result = 0.0
        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        qmk = q - k
        
        ek, ubov = GLSW.eigensystem(k)
        eqmk, ubov = GLSW.eigensystem(qmk)
        
        for band2 in range(2*num_sub):
            for band3 in range(2*num_sub):
                        
                vd = vertex.V2_cubic(band1, band2, band3, q, k, qmk)
                tmp = vd.conj() * vd/(omega - ek[band2] - eqmk[band3] + 1j*cgf)
                result += tmp
        
        ff[0] = np.real(result)
        ff[1] = np.imag(result)
        
        return 0
    
    print_header('Vegas')
    res = pycuba.Vegas(integrand_fun1, NDIM, verbose=2, ncomp=2)
    print_results('Vegas', res)
    Res = res.get('results')[0]
    integrals = Res.get('integral')
    
    return integrals
    
    
        
        
            
        
                
                
                
        
        
