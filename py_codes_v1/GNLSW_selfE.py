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

cgf = cf.convergence   # the convergence factor for the self-energy
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
    

# cuba variables

NDIM = 3
EPSREAL = 1e-3
EPSABS = 1e-3
MINEVEL = 0
MAXEVAL = 50000


def Sigma_decay(q, omega, band1):
    """
    calculates the self energy from the decay channel
    """
    
    # define the integrand as a nested function
    def integrand_decay(ndim, xx, ncopm, ff, userdata):
        """
        the decay integrand
        """
        
        result = 0.0   # the value of the integrand
        
        # integration variable: k input variable: q
        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        qmk = q - k
        
        eq, ubov_q = GLSW.eigensystem(q)
        ek, ubov_k = GLSW.eigensystem(k)
        eqmk, ubov_qmk = GLSW.eigensystem(qmk)
        
        tmp, ubov_mq = GLSW.eigensystem(-q)
        tmp, ubov_mk = GLSW.eigensystem(-k)
        tmp, ubov_mqmk = GLSW.eigensystem(-qmk)
        
        for band2 in range(2*num_sub):
            for band3 in range(2*num_sub):
                        
                vd = vertex.V2_cubic(band1, band2, band3, q, k, qmk, \
                                     ubov_q, ubov_k, ubov_qmk, \
                                     ubov_mq, ubov_mk, ubov_mqmk)
                
                tmp = vd.conj() * vd/(omega - ek[band2] - eqmk[band3] + 1j*cgf)
                result += tmp
        
        ff[0] = np.real(result)
        ff[1] = np.imag(result)
        
        return 0
    
    print_header('Vegas')

    """
     vegas:
     pycuba.Vegas(integrand, ndim, userdata, epsrel=0.001, epsabs=1e-12, \
                 verbose=0, ncomp=1, seed=None, minevel=0, maxeval=50000, \
                 nstart=1000, nincrease=500, nbatch=1000, gridno=0, statefile=, \
                 nvec=1)
    """

    #res = pycuba.Vegas(integrand_decay, NDIM, epsabs = EPSABS, \
    #                  verbose=2, ncomp=2, maxeval=1000)
    
    """
    
    The return value is always a dictionary:
 { 'neval': number of evaluations,
  'fail': 0 or 1,
  'comp': number of components, usually 1,
  'nregions': number of regions used,
  'results': [ # a list of results for each component
     {
       'integral': value of the integral,
       'error':  uncertainty,
       'prob': probability (see manual),
     },
     ...
  ]
 }
  """

    print_results('Vegas', pycuba.Vegas(integrand_decay, NDIM, epsabs = EPSABS, \
                           verbose=2, ncomp=2, maxeval=1000))
    # get values of the keyword 'result' (in the dictionary)
    # first element is the result of the integration
   # rres = res.get('results')[0] 
    #iintegral = rres.get('integral')
    iintegral = 0.0
    
    return iintegral
    
    
        
        
            
        
                
                
                
        
        
