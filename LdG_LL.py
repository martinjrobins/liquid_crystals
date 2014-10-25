# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

from dolfin import *
import numpy as np

# Create mesh

rot_step = 2*pi/20
T = 0.05
N = 50
N_b = 10**4
out_dir = 'out/LdG_LL'

mesh = UnitSquare(N,N)
def bottom(x, on_boundary):
    return near(x[1],0.0) and on_boundary
def top(x, on_boundary):
    return near(x[1],1.0) and on_boundary
def left(x, on_boundary):
    return near(x[0],0.0) and on_boundary
def right(x, on_boundary):
    return near(x[0],1.0) and on_boundary

# Define function space
V = FunctionSpace(mesh, "CG", 1)

# Define boundary condition
theta_bottom = Constant(0)
theta_top = Constant(0)
theta_left = Constant(pi/2)
theta_right = Constant(pi/2)
bc_bottom = DirichletBC(V, theta_bottom, bottom)
bc_top = DirichletBC(V, theta_top, top)
bc_left = DirichletBC(V, theta_left, left)
bc_right = DirichletBC(V, theta_right, right)
bc = [bc_bottom,bc_top,bc_left,bc_right]


theta = TrialFunction(V)
v = TestFunction(V)

# Define variational problem
a = inner(grad(theta), grad(v))*dx
L = Constant(0)*v*dx

# Compute solution
theta = Function(V)
solve(a == L, theta, bc)


# Define nonlinear problem
def trap(t,d):
    if t < d and t >= 0:
        return t/d
    elif t < 1-d:
        return 1.0
    elif t <= 1:
        return (1-t)/d
    else:
        return 0

eps = 0.02
class sExpression(Expression):
    def eval(self, value, x):
        if near(x[0],1.0) or near(x[0],0.0):
            value[0] = trap(x[1],3*eps)
        elif near(x[1],1.0) or near(x[1],0.0):
            value[0] = trap(x[0],3*eps)
        else:
            value[0] = 1.0
        
s = sExpression()

# Define function space
V = VectorFunctionSpace(mesh, "CG", 1)

# Define boundary condition
u_bottom = as_vector([cos(2*theta_bottom),sin(2*theta_bottom)])*s
u_top = as_vector([cos(2*theta_top),sin(2*theta_top)])*s
u_left = as_vector([cos(2*theta_left),sin(2*theta_left)])*s
u_right = as_vector([cos(2*theta_right),sin(2*theta_right)])*s

bc_bottom = DirichletBC(V, u_bottom, bottom)
bc_top = DirichletBC(V, u_top, top)
bc_left = DirichletBC(V, u_left, left)
bc_right = DirichletBC(V, u_right, right)
bc = [bc_bottom,bc_top,bc_left,bc_right]

Q = project(as_vector([cos(2*theta),sin(2*theta)]),V=V,bcs = bc)
file = File("%s/LdG_init.pvd"%out_dir)
file << Q
v = TestFunction(V)

#F = inner(grad(Q),grad(v))*dx() + Constant(2/eps**2)*(dot(Q,Q) - Constant(1))*elem_mult(Q,v)*dx()
F = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
    inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()

# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})


file = File("%s/LdG_solution.pvd"%out_dir)
file << Q

c = 3*eps
psep = 1.0/N
overlap = 2*psep

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

index_to_particle = np.empty([N+1,N+1],dtype=int)
particles = Particles()
count = 0
for i in range(N+1):
    for j in range(N+1):
        p = Particle()
        x = [i*psep,j*psep,0]
        p.position = Vect3d(x)
        index_to_particle[i][j] = -1
        #if not in_corners_plus_overlap(x):
        #    continue
        if in_corners(x):
            p.fixed = False
        else:
            p.fixed = True
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
        p.orientation = Vect3d(cos(p.theta),sin(p.theta),0)
        particles.append(p)
        index_to_particle[i][j] = count
        count = count + 1
     
    
U = LabwohlLasherPotential(epsilon=1,lattice_spacing=1)

list_of_overlap_particles = []
list_of_overlap_vertices = []

for i in range(mesh.num_vertices()):
    x = mesh.coordinates()[i]
    lattice_i = int(x[0]/psep + 0.5)
    lattice_j = int(x[1]/psep + 0.5)
    particles_i = index_to_particle[lattice_i][lattice_j]
    if particles_i == -1:
        continue
    #if in_corners_plus_overlap(x) and not in_corners(x):
    #if in_corners_plus_overlap(x):
    list_of_overlap_particles.append([particles_i, i])
    if in_corners(x):
        list_of_overlap_vertices.append([particles_i, i])

    
#set lattice particles according to continuum verticies
Q11_Function = Q.sub(0,deepcopy=True).vector().array()
Q12_Function = Q.sub(1,deepcopy=True).vector().array()
print 51*51,size(Q11_Function),size(Q.vector().array())

for i in list_of_overlap_particles:
    particle_i = i[0]
    vertex_i = i[1]
    print particles[particle_i].position,mesh.coordinates()[vertex_i]
    print 
#    Q11 = Q11_Function.vector()[vertex_i]
#    Q12 = Q12_Function.vector()[vertex_i]
    Q11 = Q.vector()[2*vertex_i]
    Q12 = Q.vector()[2*vertex_i+1]
    s = sqrt(Q11**2+Q12**2)
    theta = acos(Q11/s)/2.0
    #print 'theta = ',acos(Q11/s)/2.0
    particles[particle_i].theta = theta
    particles[particle_i].orientation = Vect3d(cos(theta),sin(theta),0)

    
v = particles.get_grid()
w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/LL_init.vtu'%(out_dir))
w.write()
            
#run lattice monte carlo
#f = open('%s/U.dat'%(out_dir), 'w')
#for batch in range(5):
#    tau = monte_carlo_timestep(N_b,particles,U,params)
#    print tau
#    f.write('%d %f\n'%(batch,tau))
#    f.flush()
#    v = particles.get_grid()
#    w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/LL_batch%04d.vtu'%(out_dir,batch))    
#    w.write()
    
    




