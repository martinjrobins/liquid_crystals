# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

from dolfin import *

# Create mesh
mesh = UnitSquareMesh(10, 10)
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
u_bottom = Constant(0)
u_top = Constant(0)
u_left = Constant(pi/2)
u_right = Constant(pi/2)
bc_bottom = DirichletBC(V, u_bottom, bottom)
bc_top = DirichletBC(V, u_top, top)
bc_left = DirichletBC(V, u_left, left)
bc_right = DirichletBC(V, u_right, right)
bc = [bc_bottom,bc_top,bc_left,bc_right]


theta = TrialFunction(V)
v = TestFunction(V)

# Define variational problem
a = inner(grad(theta), grad(v))*dx
L = Constant(0)*v*dx

# Compute solution
theta = Function(V)
solve(a == L, theta, bc)

# Plot solution
plot(theta, interactive=True)


# Define variational problem
#(Q11,Q12) = TrialFunction(V)
#(v11,v12) = TestFunction(V)
#
#eps = Constant(1.0)
#a = (inner(grad(Q11), grad(v11))*dx() + (2/eps**2)*(Q11*Q11 + Q12*Q12 - 1)*Q11*v11)*dx()
