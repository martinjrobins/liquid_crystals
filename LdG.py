# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

from dolfin import *

# Create mesh
mesh = UnitSquare(100, 100)
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

# Plot solution
plot(theta, interactive=True)

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


class gExpression(Expression):
    def eval(self, value, x):
        if near(x[0],1.0) or near(x[0],0.0):
            value[0] = -trap(x[1],0.06)
            value[1] = 0.0
        elif near(x[1],1.0) or near(x[1],0.0):
            value[0] = trap(x[0],0.06)
            value[1] = 0.0
        else:
            value[0] = 1.0
            value[1] = 0.0
    def value_shape(self):
        return (2,)
        
g = gExpression()
s = sqrt(dot(g,g))
# Define function space
V = VectorFunctionSpace(mesh, "CG", 1)

# Define boundary condition
u_bottom = as_vector([cos(2*theta_bottom),sin(2*theta_bottom)])
u_top = as_vector([cos(2*theta_top),sin(2*theta_top)])
u_left = as_vector([cos(2*theta_left),sin(2*theta_left)])
u_right = as_vector([cos(2*theta_right),sin(2*theta_right)])
bc_bottom = DirichletBC(V, u_bottom, bottom)
bc_top = DirichletBC(V, u_top, top)
bc_left = DirichletBC(V, u_left, left)
bc_right = DirichletBC(V, u_right, right)
bc = [bc_bottom,bc_top,bc_left,bc_right]

Q = project(as_vector([cos(2*theta),sin(2*theta)]),V=V,bcs = bc)
plot(Q)
interactive()
v = TestFunction(V)
eps = 0.02
#F = inner(grad(Q),grad(v))*dx() + Constant(2/eps**2)*(dot(Q,Q) - Constant(1))*elem_mult(Q,v)*dx()
F = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
    inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()




# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
(Q11, Q12) = Q.split()

# Plot sigma and u
plot(Q11,title='Q11')
plot(Q12,title='Q12')

s = sqrt(0.5*(Q11*Q11 + Q12*Q12))
plot(s,title='s')
theta = acos(Q11)/2.0
#sin_theta = project(Q12/(2*s*cos_theta))
#plot(as_vector([cos_theta,sin_theta]),title='director')
plot(as_vector([cos(theta),sin(theta)]))
interactive()

