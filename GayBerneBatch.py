
# In[1]:

from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi,cos,sin
from png_show import png_show


# In[2]:

L = 50.0
k = 3.0
sigma_s = 0.5
rho = 0.3
rot_step = 2*pi/50
diff_step = sigma_s/50
T = 3.2
area = (1.0/4.0)*pi*k*sigma_s**2
N = int(rho*L**2/area)

hgo_steps = 10**4
N_b = 10**4
tau_s = 10**(-4)

params = Params()
params['Dtrans'] = diff_step*50
params['Drot'] = rot_step*50
params['Temp'] = T

params['L'] = L
out_dir = 'out/GB'



# In[3]:

particles = Particles()
for i in range(N):
    p = Particle()
    eps = sigma_s
    p.position = Vect3d(uniform(eps,L-eps),uniform(eps,L-eps),0)
    theta = uniform(0,2*pi)
    u = [cos(theta),sin(theta),0]
    p.orientation = Vect3d(u)
    p.fixed = False
    particles.append(p)
    
sigma_e = sigma_s*k
n_side = int(L/sigma_e)
spacing = L/n_side
for i in range(n_side):
    p = Particle()
    p.position = Vect3d(0,(i+0.5)*spacing,0)
    p.orientation = Vect3d(0,1,0)
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d(L,(i+0.5)*spacing,0)
    p.orientation = Vect3d(0,1,0)
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d((i+0.5)*spacing,0,0)
    p.orientation = Vect3d(1,0,0)
    p.fixed = True
    particles.append(p)
    
    p = Particle()
    p.position = Vect3d((i+0.5)*spacing,L,0)
    p.orientation = Vect3d(1,0,0)
    p.fixed = True
    particles.append(p)
 
U_hgo = HGOPotential(sigma_s=sigma_s,k=k)
U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)



v = particles.get_grid()
v._get_point_data().set_active_normals('orientation')
ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
domain_source = tvtk.CubeSource(x_length=L,y_length=L,z_length=0,center=[L/2.0,L/2.0,0])
glyph = tvtk.Glyph3D(source=ellipse_source.output,input=v,scale_factor=1,vector_mode=True,orient=True)

m = tvtk.PolyDataMapper(input=glyph.output)
m2 = tvtk.PolyDataMapper(input=domain_source.output)
p = tvtk.Property(representation='w')
a2 = tvtk.Actor(mapper=m2,property=p)
a = tvtk.Actor(mapper=m)
ren = tvtk.Renderer()
ren.add_actor(a)
ren.add_actor(a2)
ren.reset_camera(L*0.1,L*0.9,L*0.1,L*0.9,-0.1,0.1)


tau = monte_carlo_timestep(hgo_steps,particles,U_hgo,params)
v = particles.get_grid()
glyph.modified()

run = 0
png_show(ren,filename='%s/init%04d'%(out_dir,run),width=800,height=800)	
f = open('%s/U%04d'%(out_dir,run), 'w')


params['Dtrans'] = diff_step
params['Drot'] = rot_step

for batch in range(50):
    tau = monte_carlo_timestep(N_b,particles,U,params)
    print tau
    f.write('%d %f\n'%(batch,tau))
    f.flush()
    v = particles.get_grid()
    glyph.modified()
    png_show(ren,filename='%s/batch%04d'%(out_dir,batch),width=800,height=800)
    
f.close()
  