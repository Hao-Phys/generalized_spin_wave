#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:42:02 2019

@author: Hao
"""

import numpy as np
import pycuba

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

def user_dat(aa):
    
      def Integrand(ndim, xx, ncomp, ff, userdata):
          x,y,z = [xx[i] for i in range(ndim.contents.value)]
          result = aa*np.sin(x)*np.cos(y)*np.exp(z)
          ff[0] = result
          return 0
      
      print_results('Vegas', pycuba.Vegas(Integrand, 3, verbose=2))
      
user_dat(12.0)
      
        