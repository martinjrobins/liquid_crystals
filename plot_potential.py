# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 11:10:54 2014

@author: robinsonm
"""

from particleSimulation import *
import numpy as np
import pylab as pl


k = 3.0
sigma_s = 0.5

 
U_hgo = HGOPotential(sigma_s=sigma_s,k=k)
U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)

r = np.arange(sigma_s/100,3*sigma_s*k,sigma_s/100)
r2 = [Vect3d(rx,0,0) for rx in r]
r1 = [Vect3d(0,0,0) for rx in r]
u1 = [Vect3d(0,1,0) for rx in r]
u2 = [Vect3d(1,0,0) for rx in r]

p = np.array([U.evaluate(rr1,uu1,rr2,uu2) for (rr1,uu1,rr2,uu2) in zip(r1,u1,r2,u2)])
phgo = np.array([U_hgo.evaluate(rr1,uu1,rr2,uu2) for (rr1,uu1,rr2,uu2) in zip(r1,u1,r2,u2)])

pl.plot(r,p,'b')
pl.hold(True)
pl.plot(r,phgo,'r')
pl.ylim(-5,5)
pl.show()