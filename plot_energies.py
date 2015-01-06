# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 02:15:41 2015

@author: mrobins
"""

import matplotlib.pyplot as plt
import numpy as np
import sys


assert len(sys.argv)==2
out_dir = sys.argv[1]

plt.figure(figsize=(6,4.5))
plt.xlabel('D ($\mu m$)')
plt.ylabel('E')

data = np.loadtxt('%s/energies.txt'%out_dir)

plt.scatter(data[:,0],data[:,1],label='D1',c='b')
plt.scatter(data[:,0],data[:,2],label='D2',c='b')
plt.scatter(data[:,0],data[:,3],label='R1',c='r')
plt.scatter(data[:,0],data[:,4],label='R2',c='r')
plt.scatter(data[:,0],data[:,5],label='R3',c='r')
plt.scatter(data[:,0],data[:,6],label='R4',c='r')
plt.legend(loc='upper left')

plt.savefig('%s/energies.pdf'%out_dir)



