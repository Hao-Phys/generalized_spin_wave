#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 19:27:33 2019

@author: hao
"""
import numpy as np

def fun_sm(alpha1, alpha2, theta, phi):
    value = np.sqrt(2)*np.exp(1j*alpha1) \
            * np.sin(theta)*np.cos(theta)*np.cos(phi) \
            + np.sqrt(2)*np.exp(-1j*alpha2) \
            * np.sin(theta)*np.cos(theta)*np.sin(phi)
    return value

def fun_smsm(alpha1, alpha2, theta, phi):
    value = np.exp(1j*(alpha1-alpha2))*(np.sin(theta))**2*np.sin(2*phi)
    return value

def fun_smsz(alpha1, alpha2, theta, phi):
    value = np.sqrt(2) \
            * (np.exp(1j*alpha1) \
            * np.sin(theta)*np.cos(theta)*np.cos(phi) \
            - np.exp(-1j*alpha2) \
            * np.sin(theta)*np.cos(theta)*np.sin(phi))
    return value

def fun_sp(alpha1, alpha2, theta, phi):
    value = np.sqrt(2)*np.exp(1j*alpha2) \
            * np.sin(theta)*np.cos(theta)*np.sin(phi) \
            + np.sqrt(2)*np.exp(-1j*alpha1) \
            * np.sin(theta)*np.cos(theta)*np.cos(phi)
    return value

def fun_spsp(alpha1, alpha2, theta, phi):
    value = np.exp(-1j*(alpha1-alpha2))*(np.sin(theta))**2*np.sin(2*phi)
    return value

def fun_spsz(alpha1, alpha2, theta, phi):
    value = np.sqrt(2) \
            *(np.exp(-1j*alpha1)*np.sin(theta)*np.cos(theta)*np.cos(phi) \
            - np.exp(1j*alpha2)*np.sin(theta)*np.cos(theta)*np.sin(phi))
    return value

def fun_sz(theta, phi):
    value = (np.sin(theta))**2*np.cos(2*phi)
    return value

def fun_sz_sq(theta):
    value = (np.sin(theta))**2
    return value


    