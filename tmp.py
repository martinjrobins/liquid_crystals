# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 00:48:17 2014

@author: mrobins
"""

import glob
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from random import sample
from utilities import import_columns

out_dir = 'out/LL/temp_calc5'
files = glob.glob('%s/U_*'%out_dir)
files.sort()
plt.figure(figsize=(6,4.5))
plt.xlabel('x')
plt.ylabel('y')
order = []
for filename in files:
    print 'doing ',filename
    data = np.loadtxt(filename)
    order.append(np.mean(data[:,0]))
    #with np.load(filename) as data:

plt.plot(range(len(files)),order)
plt.show()
    
