#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 13:42:04 2022

yesy ffind_peaks

@author: boeglinw
"""

import numpy as np
import ffind_peaks as FP

import time

# examples for testing
n_p = 10000000

x = np.linspace(0., n_p*2*np.pi, 20*n_p)
series = np.sin(x)


results = np.zeros((2,), dtype = 'int32')

pmin = np.zeros(int(len(series)/5, ), dtype='int32')
pmax = np.zeros(int(len(series)/5, ), dtype='int32')
#%%
start = time.time()

FP.find_peaks(len(series), 0.2, series, results, pmin, pmax)

end = time.time()

print(f'It took {end - start} s')