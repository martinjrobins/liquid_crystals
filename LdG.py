# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 14:17:38 2014

@author: robinsonm
"""

from dolfin import *
from setupLdG import setupLdG

eps = 0.02
mesh = UnitSquare(50,50)

#D1,D2,R1,R2,R3,R4
arg1 = {pi/2, pi/2, pi/2, pi/2, 3*pi/2, pi/2}
arg2 = {pi/2, pi/2, pi/2, pi/2, pi/2, 3*pi/2}
arg3 = {0, pi, pi, 0, pi, pi}
arg4 = {0, pi, 0, pi, pi, pi}

#D1
(F,bc,Q) = setupLdG(mesh,pi/2,pi/2,0,0,eps)



out_dir = 'out/LdG'


# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
#(Q11, Q12) = Q.split()

# Plot sigma and u
#plot(Q11,title='Q11')
#plot(Q12,title='Q12')

#s = sqrt(0.5*(Q11*Q11 + Q12*Q12))
#plot(s,title='s')
#theta = acos(Q11)/2.0
#sin_theta = project(Q12/(2*s*cos_theta))
#plot(as_vector([cos_theta,sin_theta]),title='director')
#plot(as_vector([cos(theta),sin(theta)]))
#interactive()

file = File("%s/LdG_solution_D1.pvd"%out_dir)
file << Q

(F,bc,Q) = setupLdG(mesh,pi/2,pi/2,pi,0,eps)

# Compute solution
solve(F == 0, Q, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})


file = File("%s/LdG_solution_R1.pvd"%out_dir)
file << Q
