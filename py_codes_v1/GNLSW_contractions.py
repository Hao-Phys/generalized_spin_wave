#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : GNLSW_contractions.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.13.2020
# Last Modified Date: 01.14.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import GLSW
import pycuba

num_bond = cf.num_bond
num_sub = cf.num_sub
delta_ij = cf.delta_ij
sub_idx = cf.sub_idx

NDIM = 3  # dimension of integration

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



def contractions(bond):

    bond_vec = delta_ij[:, bond]
    sub1 = sub_idx[bond, 0]
    sub2 = sub_idx[bond, 1]

    def integrand_Nij(ndim, xx, ncomp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        cphase = np.exp(-1j*2.0*np.pi*k @ bond_vec)

        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov21_c = ubov21.conj()

        counter = 0

        for flavor1 in range(2):
            for flavor2 in range(2):

                X1 = num_sub*flavor1 + sub1
                X2 = num_sub*flavor2 + sub2

                Nii = (ubov11[X1, :] * ubov11_c[X1, :]).sum()
                Njj = (ubov11[X2, :] * ubov11_c[X2, :]).sum()
                Nij = (ubov11[X1, :] * ubov11_c[X2, :] * cphase).sum()

                ff[6*counter]   = np.real(Nii)
                ff[6*counter+1] = np.imag(Nii)
                ff[6*counter+2] = np.real(Njj)
                ff[6*counter+3] = np.imag(Njj)
                ff[6*counter+4] = np.real(Nij)
                ff[6*counter+5] = np.imag(Nij)
                
                counter += 1

        return 0


    print_header('Vegas')
    
    res = pycuba.Vegas(integrand_Nij, NDIM, epsabs=1e-4, \
                       verbose=2, ncomp=24, maxeval=10000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')
    
    return iintegrals


for bond in range(1):
    res1 = contractions(bond)
    print(res1)



