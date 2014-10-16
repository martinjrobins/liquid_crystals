# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

from dolfin import *

# Create mesh
mesh = UnitSquareMesh(8, 8)
def bottom(x, on_boundary):
    return near(x[1],0) and on_boundary
def top(x, on_boundary):
    return near(x[1],L) and on_boundary
def left(x, on_boundary):
    return near(x[0],0) and on_boundary
def right(x, on_boundary):
    return near(x[0],L) and on_boundary

# Define function space
V = FunctionSpace(mesh, "Lagrange", 1)

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

# Define variational problem
#(Q11,Q12) = TrialFunction(V)
#(v11,v12) = TestFunction(V)
#
#eps = Constant(1.0)
#a = (inner(grad(Q11), grad(v11))*dx() + (2/eps**2)*(Q11*Q11 + Q12*Q12 - 1)*Q11*v11)*dx()

a2 = Constant("triangle")
c2 =  Constant("triangle")


L = f*v*dx() + g*v*ds()

# Define function for the solution
u = Function(V)

# Define goal functional (quantity of interest)
M = u*dx()

# Define error tolerance
tol = 1.e-5

# Solve equation a = L with respect to u and the given boundary
# conditions, such that the estimated error (measured in M) is less
# than tol
problem = LinearVariationalProblem(a, L, u, bc)
solver = AdaptiveLinearVariationalSolver(problem, M)
solver.parameters["error_control"]["dual_variational_solver"]["linear_solver"] = "cg"
solver.solve(tol)

solver.summary()

# Plot solution(s)
plot(u.root_node(), title="Solution on initial mesh")
plot(u.leaf_node(), title="Solution on final mesh")
interactive()