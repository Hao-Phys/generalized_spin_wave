#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : GNLSW_contractions.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.13.2020
# Last Modified Date: 01.17.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import GLSW
import pycuba

num_bond = 12
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

# contractions on site
def N_onite():

    print('---- start calculating the on-site contrations N_ii ----')

    def integrand_N_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # N_onsite is pure real, imaginary part zero
                    N_onsite = (ubov11[X1, :] * ubov11_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(N_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_N_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)

    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals
# write on-site contraction to file

# iintegrals = N_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/N_onsite.npy'
# np.save(fname, data1)


def barN_onsite():

    print('---- start calculating the on-site contrations bar_N_ii ----')

    def integrand_barN_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of barN_onsite must be real, because
                    # it corresponds to local boson number reduction
                    barN_onsite = (ubov21[X1, :] * ubov21_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(barN_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_barN_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = barN_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/barN_onsite.npy'
# np.save(fname, data1)



def Delta_onsite():

    print('---- start calculating the on-site contrations Delta_ii ----')

    def integrand_Delta_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of Delta_onsite is real
                    # need to think about the reason
                    Delta_onsite = (ubov11[X1, :] * ubov21_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(Delta_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_Delta_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = Delta_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/Delta_onsite.npy'
# np.save(fname, data1)


def bar_Delta_onsite():

    print('---- start calculating the on-site contrations Delta_ii ----')

    def integrand_barDelta_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of Delta_onsite is real
                    # need to think about the reason
                    barDelta_onsite = (ubov21[X1, :] * ubov11_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(barDelta_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_Delta_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = barDelta_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/barDelta_onsite.npy'
# np.save(fname, data1)


def N_bond():

    print('---- start calculating the bond contrations N_ij')

    def integrand_N(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for bond in range(num_bond):

            counter = 0

            bond_vec = delta_ij[:, bond]
            cphase = np.exp(-1j*2.0*np.pi*k @ bond_vec)
            sub1 = sub_idx[bond, 0]
            sub2 = sub_idx[bond, 1]

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub1
                    X2 = num_sub*flavor2 + sub2

                    # N_onsite is pure real, imaginary part zero
                    N_ij = (ubov11[X1, :] * ubov11_c[X2, :] * cphase).sum()
                    ff[8*bond+2*counter] = np.real(N_ij)
                    ff[8*bond+2*counter+1] = np.imag(N_ij)
                    counter += 1 

        return 0

    NCOMP = num_bond*8

    print_header('Vegas')

    res = pycuba.Vegas(integrand_N, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=NCOMP, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)

    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals
# write on-site contraction to file

iintegrals = N_bond()
print(iintegrals)
data1 = np.reshape(iintegrals, (num_bond, 2, 4))
fname = 'res_contractions/N_bond.npy'
np.save(fname, data1)


def barN_bond():

    print('---- start calculating the on-site contrations bar_N_ii ----')

    def integrand_barN_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of barN_onsite must be real, because
                    # it corresponds to local boson number reduction
                    barN_onsite = (ubov21[X1, :] * ubov21_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(barN_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_barN_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = barN_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/barN_onsite.npy'
# np.save(fname, data1)



def Delta_onsite():

    print('---- start calculating the on-site contrations Delta_ii ----')

    def integrand_Delta_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of Delta_onsite is real
                    # need to think about the reason
                    Delta_onsite = (ubov11[X1, :] * ubov21_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(Delta_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_Delta_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = Delta_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/Delta_onsite.npy'
# np.save(fname, data1)


def bar_Delta_onsite():

    print('---- start calculating the on-site contrations Delta_ii ----')

    def integrand_barDelta_onsite(ndim, xx, comp, ff, userdata):

        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        tmp, ubov = GLSW.eigensystem(k)
        ubov11 = ubov[:2*num_sub, :2*num_sub]
        ubov21 = ubov[2*num_sub:, :2*num_sub]
        ubov11_c = ubov11.conj()
        ubov21_c = ubov21.conj()


        for sub_lat in range(num_sub):
            counter = 0

            for flavor1 in range(2):
                for flavor2 in range(2):

                    X1 = num_sub*flavor1 + sub_lat
                    X2 = num_sub*flavor2 + sub_lat

                    # integration of Delta_onsite is real
                    # need to think about the reason
                    barDelta_onsite = (ubov21[X1, :] * ubov11_c[X2, :]).sum()
                    ff[4*sub_lat+counter] = np.real(barDelta_onsite)
                    counter += 1 

        return 0


    print_header('Vegas')

    res = pycuba.Vegas(integrand_Delta_onsite, NDIM, epsabs=1e-3, \
                       verbose=2, ncomp=16, maxeval=5000)

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')

    return iintegrals

# write on-site contraction to file
# iintegrals = barDelta_onsite()
# print(iintegrals)
# data1 = np.reshape(iintegrals, (4, 2, 2))
# fname = 'res_contractions/barDelta_onsite.npy'
# np.save(fname, data1)
