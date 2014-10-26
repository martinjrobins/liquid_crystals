# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 01:18:32 2014

@author: mrobins
"""

from dolfin import *


def setupLdG(mesh,leftbc,rightbc,bottombc,topbc,eps):
    
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
    theta_bottom = Constant(bottombc)
    theta_top = Constant(topbc)
    theta_left = Constant(leftbc)
    theta_right = Constant(rightbc)
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
    
    v = TestFunction(V)
    
    #F = inner(grad(Q),grad(v))*dx() + Constant(2/eps**2)*(dot(Q,Q) - Constant(1))*elem_mult(Q,v)*dx()
    F = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
        inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()

    return (F,bc,Q)