#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 21:12:28 2020

@author: hao
"""

import numpy as np
import sys

infile = sys.argv[1]
data = np.loadtxt(infile)

idx = data[:, 0].argsort()
data_sorted = data[idx]

outfile = infile + '.sorted'
np.savetxt(outfile, data_sorted)
