# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

import numpy as np
from dolfin import *
from setupLdG import setupLdG,calc_energy,calc_eigenvalues,average_Q_and_n_and_s
#info(NonlinearVariationalSolver.default_parameters(), True)

out_dir = 'out/LdG_bifurcation'
mesh = UnitSquareMesh(100,100)


soln_name = ['D1','D2','R1','R2','R3','R4']
leftbcs = [pi/2, pi/2, pi/2, pi/2, 3*pi/2, pi/2]
rightbcs = [pi/2, pi/2, pi/2, pi/2, pi/2, 3*pi/2]
bottombcs = [0, pi, pi, 0, pi, pi]
topbcs = [0, pi, 0, pi, pi, pi]

D = np.arange(0.02,2.0,0.02)*10**(-6)
s0 = 0.6
L = 10**(-11)
c2 = 10**6
A = sqrt(2.0*(s0**2)*c2)
eps = np.sqrt(L/(A*D**2))
print eps

f = open('%s/energies.txt'%out_dir, 'w')
fe1 = open('%s/eigenvalues1.txt'%out_dir, 'w')
fe2 = open('%s/eigenvalues2.txt'%out_dir, 'w')
fs = open('%s/stable.txt'%out_dir, 'w')
fdirector1 = open('%s/director1.txt'%out_dir, 'w')
fdirector2 = open('%s/director2.txt'%out_dir, 'w')
fq1 = open('%s/q1.txt'%out_dir, 'w')
fq2 = open('%s/q2.txt'%out_dir, 'w')
forder = open('%s/order.txt'%out_dir, 'w')


f.write('#D(x10^6)')
fe1.write('#D(x10^6)')
fe2.write('#D(x10^6)')
fs.write('#D(x10^6)')
fdirector1.write('#D(x10^6)')
fdirector2.write('#D(x10^6)')
forder.write('#D(x10^6)')
fq1.write('#D(x10^6)')
fq2.write('#D(x10^6)')
for name in soln_name:
    f.write(' %s'%name)
    fe1.write(' %s'%name)
    fe2.write(' %s'%name)
    fs.write(' %s'%name)
    fdirector1.write(' %s'%name)
    fdirector2.write(' %s'%name)
    fq1.write(' %s'%name)
    fq2.write(' %s'%name)
    forder.write(' %s'%name)
f.write('\n')
fe1.write('\n')
fe2.write('\n')
fs.write('\n')
fdirector1.write('\n')
fdirector2.write('\n')
fq1.write('\n')
fq2.write('\n')
forder.write('\n')



for (theD,theEps) in zip(D,eps):
    f.write('%f'%(theD*10**6))
    fe1.write('%f'%(theD*10**6))
    fe2.write('%f'%(theD*10**6))
    fs.write('%f'%(theD*10**6))
    fdirector1.write('%f'%(theD*10**6))
    fdirector2.write('%f'%(theD*10**6))
    fq1.write('%f'%(theD*10**6))
    fq2.write('%f'%(theD*10**6))
    forder.write('%f'%(theD*10**6))
    for (name,leftbc,rightbc,bottombc,topbc) in zip(soln_name,leftbcs,rightbcs,bottombcs,topbcs):
        print theD,name,leftbc,rightbc,bottombc,topbc
        (F,bc,Q) = setupLdG(mesh,leftbc,rightbc,bottombc,topbc,theEps)
        # Compute solution
        solve(F == 0, Q, bc,
		solver_parameters={"newton_solver":{"relative_tolerance": 1e-6,
						    "error_on_nonconvergence": False,
						    "maximum_iterations": 20}})
        file = File("%s/LdG_solution_%f_%s.pvd"%(out_dir,theD*10**6,name))
        file << Q

        f.write(' %f'%(calc_energy(Q,theEps)))

        n = 5
	eigensolver = calc_eigenvalues(mesh,leftbc,rightbc,bottombc,topbc,theEps,Q,n)

        # true if all converged eigenvalues are < 0
        stable = [(eigensolver.get_eigenvalue(i)[0] < 0) for i in range(eigensolver.get_number_converged())]
        fs.write(' %d'%(all(stable)))
        fe1.write(' %f'%(eigensolver.get_eigenvalue(0)[0]))
        fe2.write(' %f'%(eigensolver.get_eigenvalue(1)[0]))
        Q,n,s = average_Q_and_n_and_s(Q,mesh)
        fdirector1.write(' %f'%n[0])
        fdirector2.write(' %f'%n[1])
        fq1.write(' %f'%Q[0])
        fq2.write(' %f'%Q[1])
        forder.write(' %f'%s)
        print eigensolver.get_eigenvalue(0)[0]
        print eigensolver.get_eigenvalue(1)[0]
        print eigensolver.get_eigenvalue(2)[0]
        print eigensolver.get_eigenvalue(3)[0]

    f.flush()
    fe1.flush()
    fe2.flush()
    fs.flush()
    fdirector1.flush()
    fdirector2.flush()
    fq1.flush()
    fq2.flush()
    forder.flush()

    f.write('\n')
    fe1.write('\n')
    fe2.write('\n')
    fs.write('\n')
    fdirector1.write('\n')
    fdirector2.write('\n')
    fq1.write('\n')
    fq2.write('\n')
    forder.write('\n')

f.close()
fe1.close()
fe2.close()
fs.close()
fdirector1.close()
fdirector2.close()
fq1.close()
fq2.close()
forder.close()




