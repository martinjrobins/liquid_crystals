# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

import numpy as np
from setupLdG import setupLdG
from dolfin import *
import os
import sys
from math import atan2


# Create mesh

rot_step = 2*pi/20
T = 0.05
N = 100
N_LdG = 50
N_b = 10**3
out_dir = 'out/LdG_LL'
eps = 0.02

mesh = UnitSquareMesh(N_LdG,N_LdG)
(F,bc,Q) = setupLdG(mesh,pi/2,pi/2,pi,0,eps)

file = File("%s/LdG_init.pvd"%out_dir)
file << Q

# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})


file = File("%s/LdG_solution.pvd"%out_dir)
file << Q

c = 4*eps
psep = 1.0/N
psep_LdG = 1.0/N_LdG
overlap = 1.0*psep
averaging_diameter = psep_LdG


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
params['h'] = averaging_diameter
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
        
print 'added ',count,' particles...'
        
lattice_particles = Particles()
output_particles = Particles()
count = 0
for i in range(N_LdG+1):
    for j in range(N_LdG+1):
        p = Particle()
        x = [i*psep_LdG,j*psep_LdG]
        p.position = Vect3d(x[0],x[1],0)        
        p.fixed = True
        if in_corners(x):
            lattice_particles.append(p)
            count = count + 1
            p.fixed = False
        else:
            if (i==0) or (i==N_LdG):
                Q11 = cos(pi-0.0001)
                Q12 = sin(pi-0.0001)
            elif (j==0) or (j==N_LdG):
                Q11 = cos(0)
                Q12 = sin(0)
            else:
                (Q11,Q12) = Q(x)                
            s = sqrt(Q11**2+Q12**2)
            n1 = sqrt(Q11/(2*s)+0.5)
            if n1==0:
                n2 = 1.0
            else:
                n2 = Q12/(2*s*n1)
            p.theta = atan2(n2,n1)
            p.orientation = Vect3d(Q11,Q12,0)
            p.averaged_orientation = Vect3d(Q11,Q12,0)
        output_particles.append(p)
            
     
    
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
    n1 = sqrt(Q11/(2*s)+0.5)
    if n1==0:
        n2 = 1.0
    else:
        n2 = Q12/(2*s*n1)
    theta = atan2(n2,n1)
    particles[i].theta = theta
    particles[i].orientation = Vect3d(n1,n2,0)
    
v = particles.get_grid()
w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/LL_init.vtu'%(out_dir))
w.write()
            
#run lattice monte carlo
f = open('%s/U.dat'%(out_dir), 'w')
N = 5
for batch in range(N):
    if (batch==N-1):
        tau = monte_carlo_timestep(N_b,N_b,particles,lattice_particles,U,params)
    else:
        tau = monte_carlo_timestep(N_b,0,particles,lattice_particles,U,params)

    s = sqrt(tau[0]**2+tau[1]**2)
    print 's = ',s,' U = ',tau[2]
    f.write('%d %f\n'%(s,tau[2]))
    f.flush()
    
    w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/LL_batch%04d.vtu'%(out_dir,batch))    
    w.write()
    w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/LL_averaged%04d.vtu'%(out_dir,batch))    
    w.write()
    
for p in lattice_particles:
    x = [p.position[0],p.position[1]]
    i = int(p.position[0]/psep_LdG+0.5)
    j = int(p.position[1]/psep_LdG+0.5)
    index = i*(N_LdG+1)+j
    output_particles[index].orientation = p.orientation
    output_particles[index].averaged_orientation = p.averaged_orientation
    output_particles[index].theta = p.theta
    

w = tvtk.XMLUnstructuredGridWriter(input=output_particles.get_grid(), file_name='%s/finalAveraged%04d.vtu'%(out_dir,batch))    
w.write()

w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/finalBatch%04d.vtu'%(out_dir,batch))    
w.write()
    
    




