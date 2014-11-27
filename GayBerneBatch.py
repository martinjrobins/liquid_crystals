
# In[1]:

from particleSimulation import *
from tvtk.api import tvtk
from random import uniform,sample
from math import sqrt,pi,cos,sin
import os
import sys
from multiprocessing import Pool
import numpy as np



# In[2]:

L = 25.0
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


assert len(sys.argv)==2
out_dir = sys.argv[1]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
def run_simulation(run):
    print 'doing run %d...'%run
    
    params = Params()
    params['Dtrans'] = diff_step*10
    params['Drot'] = rot_step*10
    params['Temp'] = T*0.000000001
    params['h'] = averaging_diameter
    params['L'] = L
    
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
    Nl = int(L)
    for i in range(Nl+1):
        for j in range(Nl+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.fixed = True
            lattice_particles.append(p)
     
    U_hgo = HGOPotential(sigma_s=sigma_s,k=k)
    U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)
    
    tau = monte_carlo_timestep(N_b,0,particles,lattice_particles,U_hgo,params)[2]
    
    
    if not os.path.exists(out_dir+'/%04d'%run):
        os.makedirs(out_dir+'/%04d'%run)
        
    f = open('%s/%04d/U'%(out_dir,run), 'w')
    
    params['Dtrans'] = diff_step
    params['Drot'] = rot_step
    params['Temp'] = T

    N_relax = 50
    N_run = 100
    
    sample_trajectories_i = sample(range(N),100)
    sample_trajectories = np.empty([len(sample_trajectories_i),2,N_run])

    for batch in range(N_relax):
        tau = monte_carlo_timestep(N_b,0,particles,lattice_particles,U,params)        
        s = sqrt(tau[0]**2+tau[1]**2)
        print 's = ',s,' U = ',tau[2]
        f.write('%d %f %f\n'%(batch*N_b,s,tau[2]))
        f.flush()
        
    for batch in range(N_run):
        tau = monte_carlo_timestep(N_b*10,N_b,particles,lattice_particles,U,params)
        
        for i in range(len(sample_trajectories_i)):
            sample_trajectories[i,0,batch] = particles[i].position[0]    
            sample_trajectories[i,1,batch] = particles[i].position[1]
        
        s = sqrt(tau[0]**2+tau[1]**2)
        print 's = ',s,' U = ',tau[2]
        f.write('%d %f %f\n'%(batch*N_b*10+N_relax*N_b,s,tau[2]))
        f.flush()
    import numpy as np

        w = tvtk.XMLUnstructuredGridWriter(input=particles.get_grid(), file_name='%s/%04d/vtkBatch%04d.vtu'%(out_dir,run,batch))
        w.write()
    
        w = tvtk.XMLUnstructuredGridWriter(input=lattice_particles.get_grid(), file_name='%s/%04d/vtkAveraged%04d.vtu'%(out_dir,run,batch))
        w.write()
    
        np.save('%s/%04d/trajectories.npy'%(out_dir,run), sample_trajectories)    
    
        
pool = Pool(processes=4)
pool.map(run_simulation, range(4))
    
      
