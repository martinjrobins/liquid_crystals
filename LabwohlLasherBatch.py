# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from particleSimulation import *
from random import uniform
from tvtk.api import tvtk
from math import sqrt,pi,cos,sin
import os
import sys
from multiprocessing import Pool


# <codecell>

L = 50.0
rot_step = 2*pi/20
diff_step = 0
T = 0.05
N = int(L)
N_b = 10**4
tau_s = 10**(-3)

averaging_diameter = 1.01

print 'adding ',N*N,' particles...'

params = Params()
params['Dtrans'] = diff_step
params['Drot'] = rot_step
params['Temp'] = T
params['h'] = averaging_diameter
params['L'] = L

assert len(sys.argv)==2
out_dir = sys.argv[1]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
        

def run_simulation(run):
    print 'doing run %d...'%run
    
    # <codecell>
    
    particles = Particles()
    for i in range(N+1):
        for j in range(N+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.theta = 0
            p.fixed = True
            if (i==0) or (i==N):
                if (j==0) or (j==N):
                    continue
                p.theta = pi/2
            elif (j==0) or (j==N):
                if (i==0) or (i==N):
                    continue
                p.theta = 0
            else:
                p.theta = uniform(0,2*pi)
                p.fixed = False
            particles.append(p)
            
    lattice_particles = Particles()
    for i in range(N+1):
        for j in range(N+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.fixed = True
            lattice_particles.append(p)
     
    
    U = LabwohlLasherPotential(epsilon=1,lattice_spacing=1)

    f = open('%s/U%04d'%(out_dir,run), 'w')    
    for batch in range(15):
        tau = monte_carlo_timestep(N_b,N_b,particles,lattice_particles,U,params)
        s = sqrt(tau[0]**2+tau[1]**2)
        print 's = ',s,' U = ',tau[2]
        f.write('%f %f\n'%(s,tau[2])
        f.flush()

    w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/finalBatch%04d.vtu'%(out_dir,run))
    w.write()
    
    w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/finalAveraged%04d.vtu'%(out_dir,run))
    w.write()    
    
    f.close()
    
pool = Pool(processes=8)              # start 4 worker processes
pool.map(run_simulation, range(100))


    
