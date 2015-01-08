# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 16:07:12 2015

@author: robinsonm
"""

import numpy as np
from dolfin import *
from setupLdG import setupLdG,calc_energy
from particleSimulation import *

import multiprocessing



out_dir = 'out/LL_bifurcation'
mesh = UnitSquareMesh(50,50)

soln_name = ['D1','R1',]
leftbcs = [pi/2,  pi/2]
rightbcs = [pi/2,  pi/2]
bottombcs = [0, pi]
topbcs = [0, 0]


D = np.arange(0.01,2,0.01)*10**(-6)
s0 = 0.6
L = 10**(-11)
eps = (2./3.)*D*L




def run_simulation(args):
    print 'doing run D = %f, name = %s...'%(args['D']*10**6,args['name'])
    
    L = 50.0
    rot_step = 2*pi/20
    diff_step = 0
    T = 0.05
    N = int(L)
    N_b = 10**4
    
    averaging_diameter = 1.01
    
    eps = 0.02
    mesh = UnitSquareMesh(50,50)
    (F,bc,Q) = setupLdG(mesh,args['leftbc'],args['rightbc'],args['bottombc'],args['topbc'],eps)
        
    solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
    
    particles = Particles()
    for i in range(N+1):
        for j in range(N+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.theta = 0
            p.fixed = True
            if ((i==0) and (j==0)) or ((i==0) and (j==N)) or ((i==N) and (j==0)) or ((i==N) and (j==N)):
                continue
            if (i==0):
                p.theta = args['leftbc']
            elif (i==N):
                p.theta = args['rightbc']
            elif (j==0):
                p.theta = args['bottombc']
            elif (j==N):
                p.theta = args['topbc']
            else:
                x = [p.position[0]/L,p.position[1]/L]
                (Q11,Q12) = Q(x)
                s = sqrt(Q11**2+Q12**2)
                n1 = sqrt(Q11/(2*s)+0.5)
                if n1==0:
                    n2 = 1.0
                else:
                    n2 = Q12/(2*s*n1)
                p.theta = atan2(n2,n1)
                p.fixed = False
            particles.append(p)
            
    lattice_particles = Particles()
    for i in range(N+1):
        for j in range(N+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.fixed = True
            lattice_particles.append(p)
     
    
    U = LabwohlLasherPotential(epsilon=args['eps'],lattice_spacing=1)

    params = Params()
    params['Dtrans'] = diff_step
    params['Drot'] = rot_step
    params['Temp'] = T
    params['h'] = averaging_diameter
    params['L'] = L

    f = open('%s/U_%f_%s'%(out_dir,args['D']*10**6,args['name']), 'w')    
    for batch in range(10):
        tau = monte_carlo_timestep(N_b,N_b,particles,lattice_particles,U,params)
        s = sqrt(tau[0]**2+tau[1]**2)
        print 's = ',s,' U = ',tau[2]
        f.write('%f %f\n'%(s,tau[2]))
        f.flush()

    w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/finalBatch_%f_%s.vtu'%(out_dir,args['D']*10**6,args['name']))
    w.write()
    
    w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/finalAveraged_%f_%s.vtu'%(out_dir,args['D']*10**6,args['name']))
    w.write()    
    
    f.close()
    
    return tau[2]


if __name__ == '__main__':
    print eps

    pool = multiprocessing.Pool(len(soln_name))

    f = open('%s/energies.txt'%out_dir, 'w')

    f.write('#D(x10^6)')
    for name in soln_name:
        f.write(' %s'%name)
    f.write('\n')

    for (theD,theEps) in zip(D,eps):
        f.write('%f'%(theD*10**6))
        args = []
        for (name,leftbc,rightbc,bottombc,topbc) in zip(soln_name,leftbcs,rightbcs,bottombcs,topbcs):
            args.append(
                {'D': theD, 'eps': theEps, 'name': name, 
                'leftbc': leftbc, 'rightbc': rightbc, 'bottombc': bottombc, 'topbc': topbc}                    
                )
           
        for result in pool.map(run_simulation, args):
            f.write(' %f'%(result))
            
        f.write('\n')
        f.flush()
        
    f.close()
