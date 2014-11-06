
# In[1]:

from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi,cos,sin
import os
import sys


# In[2]:

L = 75.0
k = 3.0
sigma_s = 0.5
rho = 0.3
rot_step = 2*pi/25
diff_step = sigma_s/20
T = 3.2
area = (1.0/4.0)*pi*k*sigma_s**2
#N = int(rho*L**2/area)
N = int(rho*(L**2)/(sigma_s**2))


averaging_diameter = 2.5

print 'adding ',N,' particles...'


N_b = 10**5
tau_s = 10**(-4)

params = Params()
params['Dtrans'] = diff_step*10
params['Drot'] = rot_step*10
params['Temp'] = T*0.000000001
params['h'] = averaging_diameter
params['L'] = L

assert len(sys.argv)==2
out_dir = sys.argv[1]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

particles = Particles()
for i in range(N):
    p = Particle()
    eps = sigma_s
    p.position = Vect3d(uniform(eps,L-eps),uniform(eps,L-eps),0)
    p.theta = uniform(0,2*pi)
    p.fixed = False
    particles.append(p)
    
sigma_e = sigma_s*k
n_side = int(L/sigma_e)
spacing = L/n_side
for i in range(n_side):
    p = Particle()
    p.position = Vect3d(0,(i+0.5)*spacing,0)
    p.theta = pi/2
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d(L,(i+0.5)*spacing,0)
    p.theta = pi/2
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d((i+0.5)*spacing,0,0)
    p.theta = 0
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d((i+0.5)*spacing,L,0)
    p.theta = 0
    p.fixed = True
    particles.append(p)
    
lattice_particles = Particles()
N = int(L)
for i in range(N+1):
    for j in range(N+1):
        p = Particle()
        p.position = Vect3d(i,j,0)
        p.fixed = True
        lattice_particles.append(p)
 
U_hgo = HGOPotential(sigma_s=sigma_s,k=k)
U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)

tau = monte_carlo_timestep(N_b,N_b/10,particles,lattice_particles,U_hgo,params)[2]

run = 0
f = open('%s/U%04d'%(out_dir,run), 'w')


params['Dtrans'] = diff_step
params['Drot'] = rot_step
params['Temp'] = T

for batch in range(200):
    tau = monte_carlo_timestep(N_b,N_b/10,particles,lattice_particles,U,params)[2]
    print tau
    f.write('%d %f\n'%(batch,tau))
    f.flush()

    w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/vtkBatch%04d.vtu'%(out_dir,batch))
    w.write()
    
    w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/vtkAveraged%04d.vtu'%(out_dir,batch))
    w.write()
    
f.close()

  
