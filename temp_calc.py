# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:32:15 2015

@author: robinsonm
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 16:07:12 2015

@author: robinsonm
"""

import numpy as np
from particleSimulation import *
import multiprocessing
from random import uniform
from tvtk.api import tvtk
from math import sqrt,pi,cos,sin
import os
import sys



out_dir = 'out/LL/temp_calc2'


def run_simulation(args):
    print 'doing T = %f...'%(args['T'])
    
    L = 50.0
    rot_step = pi
    diff_step = 0
    T = args['T']
    N = int(L)
    N_b = 10**4
    
    averaging_diameter = 1.01
    
    
    particles = Particles()
    for i in range(N):
        for j in range(N):
            p = Particle()
            p.position = Vect3d(i+0.5,j+0.5,0)
            p.theta = uniform(0,2*pi)
            p.fixed = False
            particles.append(p)
            
    lattice_particles = Particles()
    for i in range(N):
        for j in range(N):
            p = Particle()
            p.position = Vect3d(i+0.5,j+0.5,0)
            p.fixed = True
            lattice_particles.append(p)
     
    
    U = LabwohlLasherPotential(epsilon=1.,lattice_spacing=1)

    params = Params()
    params['Dtrans'] = diff_step
    params['Drot'] = rot_step
    params['Temp'] = T
    params['h'] = averaging_diameter
    params['L'] = L
    params['periodic'] = True

    f = open('%s/U_%f'%(out_dir,args['T']), 'w')    
    for batch in range(10):
        tau = monte_carlo_timestep(N_b,N_b,particles,lattice_particles,U,params)
        s = sqrt(tau[0]**2+tau[1]**2)
        print 's = ',s,' U = ',tau[2]
        f.write('%f %f\n'%(s,tau[2]))
        f.flush()

    w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/finalBatch_%f.vtu'%(out_dir,args['T']))
    w.write()
    
    w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/finalAveraged_%f.vtu'%(out_dir,args['T']))
    w.write()    
    
    f.close()
    
    return s


if __name__ == '__main__':

    cpus = 10
    pool = multiprocessing.Pool(cpus)

    f = open('%s/search_history.txt'%out_dir, 'w')

    T = [10**x for x in np.arange(-1,1,0.1)]
    args = [{'T':theT} for theT in T]
    
    for theT,result in zip(T,pool.map(run_simulation, args)):
        f.write('%f %f\n'%(theT,result))
        f.flush()
        
    f.close()
