# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

import numpy as np
from setupLdG import setupLdG
from dolfin import *


# Create mesh

rot_step = 2*pi/20
T = 0.05
N = 100
N_LdG = 50
N_b = 10**4
out_dir = 'out/LdG_LL'
eps = 0.02

mesh = UnitSquare(N_LdG,N_LdG)
(F,bc,Q) = setupLdG(mesh,pi/2,pi/2,0,0,eps)

file = File("%s/LdG_init.pvd"%out_dir)
file << Q

# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})


file = File("%s/LdG_solution.pvd"%out_dir)
file << Q

c = 4*eps
psep = 1.0/N
overlap = 1.0*psep

def in_corners(x):
    return True if x[0]**2 + x[1]**2 <= c**2 or \
                       x[0]**2 + (x[1]-1.0)**2 <= c**2 or \
                       (x[0]-1.0)**2 + x[1]**2 <= c**2 or \
                       (x[0]-1.0)**2 + (x[1]-1.0)**2 <= c**2 \
                else False

def in_corners_plus_overlap(x):
    return True if x[0]**2 + x[1]**2 <= (c + overlap)**2 or \
                       x[0]**2 + (x[1]-1.0)**2 <= (c + overlap)**2 or \
                       (x[0]-1.0)**2 + x[1]**2 <= (c + overlap)**2 or \
                       (x[0]-1.0)**2 + (x[1]-1.0)**2 <= (c + overlap)**2\
                else False
                
def in_corners_minus_overlap(x):
    return True if x[0]**2 + x[1]**2 <= (c - overlap)**2 or \
                       x[0]**2 + (x[1]-1.0)**2 <= (c - overlap)**2 or \
                       (x[0]-1.0)**2 + x[1]**2 <= (c - overlap)**2 or \
                       (x[0]-1.0)**2 + (x[1]-1.0)**2 <= (c - overlap)**2\
                else False


from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
#from math import sqrt,pi,cos,sin

params = Params()
params['Dtrans'] = 0
params['Drot'] = rot_step
params['Temp'] = T
params['L'] = 1.0

list_of_overlap_particles = []
particles = Particles()
count = 0
for i in range(N+1):
    for j in range(N+1):
        p = Particle()
        x = [i*psep,j*psep,0]
        p.position = Vect3d(x)
        p.fixed = False
        if not in_corners_plus_overlap(x):
            continue
        if (i==0) or (i==N):
            if (j==0) or (j==N):
                continue
            p.theta = pi/2
            p.fixed = True
        elif (j==0) or (j==N):
            if (i==0) or (i==N):
                continue
            p.theta = 0
            p.fixed = True
        else:
            p.theta = uniform(0,2*pi)
        if not in_corners(x):
            p.fixed = True
            list_of_overlap_particles.append(count)
        p.orientation = Vect3d(cos(p.theta),sin(p.theta),0)
        particles.append(p)
        count = count + 1
     
    
U = LabwohlLasherPotential(epsilon=1,lattice_spacing=psep)

list_of_overlap_vertices = []
for i in range(mesh.num_vertices()):
    x = mesh.coordinates()[i]
    if in_corners(x):
        list_of_overlap_vertices.append(i)

    
#set lattice particles according to continuum verticies
for i in list_of_overlap_particles:
    x = [particles[i].position[0],particles[i].position[1]]
    (Q11,Q12) = Q(x)
    s = sqrt(Q11**2+Q12**2)
    theta = acos(Q11/s)/2.0
    particles[i].theta = theta
    particles[i].orientation = Vect3d(cos(theta),sin(theta),0)

    
v = particles.get_grid()
w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/LL_init.vtu'%(out_dir))
w.write()
            
#run lattice monte carlo
f = open('%s/U.dat'%(out_dir), 'w')
for batch in range(5):
    tau = monte_carlo_timestep(N_b,particles,U,params)
    print tau
    f.write('%d %f\n'%(batch,tau))
    f.flush()
    v = particles.get_grid()
    w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/LL_batch%04d.vtu'%(out_dir,batch))    
    w.write()
    
    




