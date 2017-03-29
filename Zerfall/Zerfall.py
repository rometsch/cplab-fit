# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:30:22 2017

@author: thomas
"""

import numpy as np;

data = np.loadtxt("Agdecay.dat");

delta_t = 5;

gammas = data[:,1];
times = data[:,0]*delta_t;

