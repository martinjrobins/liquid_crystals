# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

import numpy as np
from dolfin import *
from setupLdG import setupLdG,calc_energy


out_dir = 'out/LdG_bifurcation'
mesh = UnitSquareMesh(50,50)


soln_name = ['D1','D2','R1','R2','R3','R4']
leftbcs = [pi/2, pi/2, pi/2, pi/2, 3*pi/2, pi/2]
rightbcs = [pi/2, pi/2, pi/2, pi/2, pi/2, 3*pi/2]
bottombcs = [0, pi, pi, 0, pi, pi]
topbcs = [0, pi, 0, pi, pi, pi]

D = np.arange(0.1,10,0.1)*10**(-6)
s0 = 0.6
L = 10**(-11)
c2 = 10**6
A = sqrt(2.0*(s0**2)*c2)
eps = np.sqrt(L/(A*D**2))
print eps

f = open('%s/energies.txt'%out_dir, 'w')

f.write('#D(x10^6)')
for name in soln_name:
    f.write(' %s'%name)
f.write('\n')

for (theD,theEps) in zip(D,eps):
    f.write('%f'%(theD*10**6))
    for (name,leftbc,rightbc,bottombc,topbc) in zip(soln_name,leftbcs,rightbcs,bottombcs,topbcs):
        print theD,name,leftbc,rightbc,bottombc,topbc
	(F,bc,Q) = setupLdG(mesh,leftbc,rightbc,bottombc,topbc,theEps)
        # Compute solution
        solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
                                        
        file = File("%s/LdG_solution_%f_%s.pvd"%(out_dir,theD*10**6,name))
        file << Q
        
        f.write(' %f'%(calc_energy(Q,theEps)))
        f.flush()
    f.write('\n')
    
f.close()

        
        
        
