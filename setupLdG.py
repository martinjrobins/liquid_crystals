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
                #value[0] = trap(x[1],3*eps)
                value[0] = trap(x[1],0.06)
            elif near(x[1],1.0) or near(x[1],0.0):
                #value[0] = trap(x[0],3*eps)
                value[0] = trap(x[0],0.06)
            else:
                value[0] = 1.0

    s = sExpression()

    # Define function space
    V = VectorFunctionSpace(mesh, "CG", 1)
    Q = Function(V,name='u')

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

    Q.assign(project(as_vector([cos(2*theta),sin(2*theta)]),V=V,bcs = bc))

    v = TestFunction(V)

    #F = inner(grad(Q),grad(v))*dx() + Constant(2/eps**2)*(dot(Q,Q) - Constant(1))*elem_mult(Q,v)*dx()
    F = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
        inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()


    return (F,bc,Q)

def calc_eigenvalues(mesh,leftbc,rightbc,bottombc,topbc,eps,Qin,n):
    def bottom(x, on_boundary):
        return near(x[1],0.0) and on_boundary
    def top(x, on_boundary):
        return near(x[1],1.0) and on_boundary
    def left(x, on_boundary):
        return near(x[0],0.0) and on_boundary
    def right(x, on_boundary):
        return near(x[0],1.0) and on_boundary

    # Define boundary condition
    theta_bottom = Constant(bottombc)
    theta_top = Constant(topbc)
    theta_left = Constant(leftbc)
    theta_right = Constant(rightbc)

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
                #value[0] = trap(x[1],3*eps)
                value[0] = trap(x[1],0.06)
            elif near(x[1],1.0) or near(x[1],0.0):
                #value[0] = trap(x[0],3*eps)
                value[0] = trap(x[0],0.06)
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

    def bdry(x,on_boundary):
        return on_boundary
    bc = DirichletBC(V,as_vector([Constant(0.0),Constant(0.0)]),bdry)

    v = TestFunction(V)
    q = TrialFunction(V)
    Q = Function(V)
    Q.assign(Qin)


    #F = inner(grad(Q),grad(v))*dx() + Constant(2/eps**2)*(dot(Q,Q) - Constant(1))*elem_mult(Q,v)*dx()
    F = -inner(grad(q[0]), grad(v[0]))*dx - (2/eps**2)*(q[0]*(3*Q[0]*Q[0] + Q[1]*Q[1] - 1) - 2*q[1]*Q[0]*Q[1])*v[0]*dx + \
        -inner(grad(q[1]), grad(v[1]))*dx - (2/eps**2)*(q[1]*(3*Q[1]*Q[1] + Q[0]*Q[0] - 1) - 2*q[0]*Q[1]*Q[0])*v[1]*dx

    L = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
        inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()
    #L = inner(Constant(1),v[0])*dx + inner(Constant(1),v[1])*dx


    #F = inner(grad(Q[0]), grad(v[0]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[0]*v[0]*dx() + \
     #   inner(grad(Q[1]), grad(v[1]))*dx() + (2/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)*Q[1]*v[1]*dx()

    m = inner(q,v)*dx


    # Assemble matrix
    #A = PETScMatrix()
    #M = PETScMatrix()
    #b = PETScVector()
    #assemble_system(F,L,bc,A_tensor=A,b_tensor=b)
    #assemble_system(m,L,bc,A_tensor=M,b_tensor=b)
    A,_ = assemble_system(F,L,bc)
    M,_ = assemble_system(m,L,bc)
    #assemble_system(m,L,bc,A_tensor=M,b_tensor=b)
    #assemble(F,tensor=A)
    #bc.apply(A)
    #assemble(m,tensor=M)
    bc.zero(M)
    A = down_cast(A)
    M = down_cast(M)
    #assemble(L,tensor=M)

    #dQ = TrialFunction(V)
    #a = derivative(F, Q, dQ)
    #assemble(a,tensor=A)
    #assemble_system(a,F,bc,A_tensor=A,b_tensor=b)

    # Create eigensolver
    eigensolver = SLEPcEigenSolver(A,M)
    # Specify the solution method (default is krylov-schur)
    eigensolver.parameters["solver"] = "krylov-schur"
    # Specify the part of the spectrum desired
    eigensolver.parameters["spectrum"] = "smallest magnitude"
    # Specify the problem type (this can make a big difference)
    eigensolver.parameters["problem_type"] = "gen_hermitian"
    # Use the shift-and-invert spectral transformation
    eigensolver.parameters["spectral_transform"] = "shift-and-invert"
    # Specify the shift
    eigensolver.parameters["spectral_shift"] = 1.0e-10

    #eigensolver.parameters['spectrum'] = 'smallest magnitude'
    #eigensolver.parameters['tolerance'] = 1e-10
    #eigensolver.parameters['maximum_iterations'] = 100

    # Compute all eigenvalues of A x = \lambda x
    eigensolver.solve(n)

    return eigensolver



def calc_energy(Q,eps):
    E = (inner(grad(Q[0]), grad(Q[0])) + inner(grad(Q[1]), grad(Q[1])) + (1/eps**2)*(Q[0]*Q[0] + Q[1]*Q[1] - 1)**2) * dx()
    return assemble(E)

def average_n_and_s(Q_field):
    Q = [assemble(Q_field[0]*dx),assemble(Q_field[1]*dx)]
    print Q
    s = sqrt(Q[0]**2 + Q[1]**2)
    n1 = sqrt(0.5*(Q[0]/s + 1))
    if Q[1]<1e-8:
        n2 = 0
    else:
        n2 = Q[1]*n1/(Q[0]+s)

    return [n1,n2],s

