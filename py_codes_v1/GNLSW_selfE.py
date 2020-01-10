#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:27:15 2019

@author: Hao
"""

import numpy as np
#import numpy.matlib
import cofig as cf
import GLSW
import GNLSW_vertex_new as vertex
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

# on-shell approximation

def Sigma_decay(q, eq, ubov_q, ubov_mq):
    """
    Calculates the decay self-energy using on-shell approximation

    Parameters
    ----------
    q : np.array((3, ))
        incoming momentum.
    eq : np.array((4*num_sub, ))
        incoming energy.
    ubov_q : np.array((4*num_sub, 4*num_sub))
        Bogoliubov matrix of q.
    ubov_mq : np.array((4*num_sub, 4*num_sub))
        Bogoliubov matrix of -q.

    Returns
    -------
    np.array((2*num_sub, ))
        Self-energies for 2*num_sub bands.

    """
      
    # define the integrand as a nested function
    def integrand_decay(ndim, xx, ncopm, ff, userdata):
        """
        Integrand passed to cuba.

        Parameters
        ----------
        ndim : int
            dimension of integration.
        xx : np.array((3, ))
            input momentum.
        ncopm : int
            number of components.
        ff : np.array()
            integrand.
        userdata : c-type
            I don't know how to use this for python. try to avoid it

        Returns
        -------
        int
            value not relevent unless -999, see mannual.

        """

                
        # integration variable: k
        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        qmk = q - k
               
        ek, ubov_k = GLSW.eigensystem(k)
        eqmk, ubov_qmk = GLSW.eigensystem(qmk)
        tmp, ubov_mk = GLSW.eigensystem(-k)
        tmp, ubov_mqmk = GLSW.eigensystem(-qmk)
        
        
        tmpmat1 = ek[:2*num_sub][:, None]
        tmpmat2 = eqmk[:2*num_sub][None, :]
        esum = tmpmat1 + tmpmat2
        #denomin = eq[:2*num_sub, None, None] - esum[None, :, :] + 1j*cgf

        denomin = eq[None, None, :2*num_sub] - esum[:, :, None] + 1j*cgf

        #v_decay = vertex.V_cubic_decay(q, k, qmk, ubov_q, ubov_k, ubov_qmk, \
        #                               ubov_mq, ubov_mk, ubov_mqmk)
        
        v_decay = vertex.V_cubic_decay(k, qmk, q, ubov_k, ubov_qmk, ubov_q, \
                                       ubov_mk, ubov_mqmk, ubov_mq)
        Intmat = (v_decay.conj() * v_decay)/denomin
        # 0.5 is the symmetry factor, return to the "on-shell" self-energy 
        # for all bands
        Intvec = 0.5*Intmat.sum(axis=(0, 1))  
        
        

        for band in range(2*num_sub):
            ff[2*band] = np.real(Intvec[band])
            ff[2*band+1] = np.imag(Intvec[band])
        
        return 0
    
    print_header('Vegas')

    """
     vegas:
     pycuba.Vegas(integrand, ndim, userdata, epsrel=0.001, epsabs=1e-12, \
                 verbose=0, ncomp=1, seed=None, minevel=0, maxeval=50000, \
                 nstart=1000, nincrease=500, nbatch=1000, gridno=0, statefile=, \
                 nvec=1)
    """

    res = pycuba.Vegas(integrand_decay, NDIM, epsabs = 1e-5, \
                      verbose=2, ncomp=4*num_sub, maxeval=10000)
    
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

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')
     
    return iintegrals    
    
    
def Sigma_source(q, eq, ubov_q, ubov_mq):
    """
    Calculates the source self-energy using on-shell approximation

    Parameters
    ----------
    q : np.array((3, ))
        incoming momentum.
    eq : np.array((4*num_sub, ))
        incoming energy.
    ubov_q : np.array((4*num_sub, 4*num_sub))
        Bogoliubov matrix of q.
    ubov_mq : np.array((4*num_sub, 4*num_sub))
        Bogoliubov matrix of -q.

    Returns
    -------
    np.array((2*num_sub, ))
        Self-energies for 2*num_sub bands.

    """
      
    # define the integrand as a nested function
    def integrand_source(ndim, xx, ncopm, ff, userdata):
        """
        Integrand passed to cuba.

        Parameters
        ----------
        ndim : int
            dimension of integration.
        xx : np.array((3, ))
            input momentum.
        ncopm : int
            number of components.
        ff : np.array()
            integrand.
        userdata : c-type
            I don't know how to use this for python. try to avoid it

        Returns
        -------
        int
            value not relevent unless -999, see mannual.

        """

                
        # integration variable: k
        k1, k2, k3 = [xx[i] for i in range(ndim.contents.value)]
        k = np.array([k1, k2, k3])
        qmk = q - k
               
        ek, ubov_k = GLSW.eigensystem(k)
        eqmk, ubov_qmk = GLSW.eigensystem(qmk)
        tmp, ubov_mk = GLSW.eigensystem(-k)
        tmp, ubov_mqmk = GLSW.eigensystem(-qmk)
        
        
        tmpmat1 = ek[:2*num_sub][:, None]
        tmpmat2 = eqmk[:2*num_sub][None, :]
        esum = tmpmat1 + tmpmat2
        denomin = eq[:2*num_sub, None, None] + esum[None, :, :] - 1j*cgf
        v_source = vertex.V_cubic_source(q, k, qmk, ubov_q, ubov_k, ubov_qmk, \
                                       ubov_mq, ubov_mk, ubov_mqmk)
        
        Intmat = (v_source.conj() * v_source)/denomin
        # -0.5 is the symmetry factor, return to the "on-shell" self-energy 
        # for all bands
        Intvec = -0.5*Intmat.sum(axis=(1, 2))  
        
        
        for band in range(2*num_sub):
            ff[band] = np.real(Intvec[band])
        
        return 0
    
    print_header('Vegas')

    """
     vegas:
     pycuba.Vegas(integrand, ndim, userdata, epsrel=0.001, epsabs=1e-12, \
                 verbose=0, ncomp=1, seed=None, minevel=0, maxeval=50000, \
                 nstart=1000, nincrease=500, nbatch=1000, gridno=0, statefile=, \
                 nvec=1)
    """

    res = pycuba.Vegas(integrand_source, NDIM, epsabs = 1e-5, \
                      verbose=2, ncomp=2*num_sub, maxeval=10000)
    
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

    print_results('Vegas', res)

    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')
    
    return iintegrals  


def greenfunction_gnlsw(omega, q, selfE_re, selfE_im):
    
    """
    calculates the green function of H-P bosons after 1/NS correction
    
    Parameters
    ----------
    omega   : np.array((anylength,))
              the frequency 
    q       : np.array((3,))
              input momentum
    selfE_re: np.array((2*num_sub, ))
              real part self-energy using on-shell approximation 
    
    selfE_im: np.array((2*num_sub, ))
              imag part self-energy using on-shell approximation
              
    Returns
    -------
    gf: np.array((4*num_sub, 4*num_sub)) 
    """
    
    len_omega = len(omega)
    gf = np.zeros([len_omega, 4*num_sub, 4*num_sub], dtype=complex)
    
    ek, ubov = GLSW.eigensystem(q)
    ek_m, ubov_m = GLSW.eigensystem(-q)
    
    ubov_m = ubov_m.conj()
    minus_mat = np.zeros([4*num_sub, 4*num_sub], dtype=complex)
    plus_mat = np.zeros([4*num_sub, 4*num_sub], dtype=complex)
    
    u11 = ubov[:2*num_sub, :2*num_sub]
    u21 = ubov[2*num_sub:, :2*num_sub]
    u11_m = ubov_m[:2*num_sub, :2*num_sub]
    u21_m = ubov_m[2*num_sub:, :2*num_sub]
    
    for band in range(2*num_sub):
        
        minus_mat[:2*num_sub, :2*num_sub] = np.outer(u11[:, band], u11[:, band].conj())
        minus_mat[:2*num_sub, 2*num_sub:] = np.outer(u11[:, band], u21[:, band].conj())
        minus_mat[2*num_sub:, :2*num_sub] = np.outer(u21[:, band], u11[:, band].conj())
        minus_mat[2*num_sub:, 2*num_sub:] = np.outer(u21[:, band], u21[:, band].conj())
        
        plus_mat[:2*num_sub, :2*num_sub] = np.outer(u21_m[:, band], u21_m[:, band].conj())
        plus_mat[:2*num_sub, 2*num_sub:] = np.outer(u21_m[:, band], u11_m[:, band].conj())
        plus_mat[2*num_sub:, :2*num_sub] = np.outer(u11_m[:, band], u21_m[:, band].conj())
        plus_mat[2*num_sub:, 2*num_sub:] = np.outer(u11_m[:, band], u11_m[:, band].conj())
        
        tmp1 = np.reshape(omega - (ek[band] + selfE_re[band]) - 1j*selfE_im[band] + 1j*0.06, len_omega)
        tmp2 = np.reshape(omega + (ek_m[band] + selfE_re[band]) - 1j*selfE_im[band] + 1j*0.06, len_omega)
        
        # use this to avoid loop 
        gf += - minus_mat/(tmp1[:, None, None]) + plus_mat/(tmp2[:, None, None])
        
    return gf



def intensity(omega, qx, qy, qz, selfE_re, selfE_im):
    
    """ 
    output neutron intensity up to 1/NS order for given qx, qy, and qz
    WARNING: ignores H-F corrections (quartic), because they are small!
    """
    
    qq = np.array([qx, qy, qz])
    q1, q2, q3 = cf.kxyTok12(qx, qy, qz)
    len_omega = len(omega)
    q = np.array([q1, q2, q3])
    
    chi_mat = np.zeros((len_omega, 3, 3), dtype=complex)
    
    gf = greenfunction_gnlsw(omega, q, selfE_re, selfE_im)
    gf11 = gf[:, :2*num_sub, :2*num_sub]
    gf12 = gf[:, :2*num_sub, 2*num_sub:]
    gf21 = gf[:, 2*num_sub:, :2*num_sub]
    gf22 = gf[:, 2*num_sub:, 2*num_sub:] 
    
    for sub1 in range(num_sub):
        for sub2 in range(num_sub):
                        
            for m in range(2):
                for mp in range(2):
                    
                    mat1, mat2, mat3, mat4 = GLSW.int_mat(sub1, sub2, m, mp)
                    
                    element1 = gf12[:, 4*m+sub1, 4*mp+sub2]
                    element2 = gf21[:, 4*m+sub1, 4*mp+sub2]
                    element3 = gf11[:, 4*m+sub1, 4*mp+sub2]
                    element4 = gf22[:, 4*m+sub1, 4*mp+sub2]
                    
                    chi_mat += (-1/num_sub)\
                            * (element1[:, None, None]*mat1 + element2[:, None, None]*mat2 \
                             + element3[:, None, None]*mat3 + element4[:, None, None]*mat4) 
                             
    sqw_mat = -2.0 * np.imag(chi_mat)
    ff = cf.formfactor(qq)
    sc_inten = ff*cf.projector(qx, qy, qz, sqw_mat)
    
    return sc_inten

    
        
            
        
                
                
                
        
        
