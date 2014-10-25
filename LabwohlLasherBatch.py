# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi,cos,sin
from png_show import png_show

# <codecell>

L = 50.0
rot_step = 2*pi/20
diff_step = 0
T = 0.05
N = 50
N_b = 10**4
tau_s = 10**(-3)

params = Params()
params['Dtrans'] = diff_step
params['Drot'] = rot_step
params['Temp'] = T
params['L'] = L
out_dir = 'out/LL'

for run in range(100):
    print 'doing run %d...'%run
    
    # <codecell>
    
    particles = Particles()
    for i in range(N+1):
        for j in range(N+1):
            p = Particle()
            p.position = Vect3d(i,j,0)
            p.theta = 0
            p.fixed = True
            if (i==0) or (i==N):
                if (j==0) or (j==N):
                    continue
                p.theta = pi/2
            elif (j==0) or (j==N):
                if (i==0) or (i==N):
                    continue
                p.theta = 0
            else:
                p.theta = uniform(0,2*pi)
                p.fixed = False
            particles.append(p)
     
    
    U = LabwohlLasherPotential(epsilon=1,lattice_spacing=1)
    
    # <codecell>
    
    v = particles.get_grid()
    v._get_point_data().set_active_normals('averaged_orientation')
    ellipse = tvtk.ParametricEllipsoid(x_radius=1.4/2.0,y_radius=0.5/2.0,z_radius=0.5/4.0)
    ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
    glyph = tvtk.Glyph3D(source=ellipse_source.output,input=v,scale_factor=1,vector_mode=True,orient=True)
    domain_source = tvtk.CubeSource(x_length=L,y_length=L,z_length=0,center=[L/2.0,L/2.0,0])
    m = tvtk.PolyDataMapper(input=glyph.output)
    a = tvtk.Actor(mapper=m)
    m2 = tvtk.PolyDataMapper(input=domain_source.output)
    p = tvtk.Property(representation='w')
    a2 = tvtk.Actor(mapper=m2,property=p)
    ren = tvtk.Renderer()
    ren.add_actor(a)
    ren.add_actor(a2)
    ren.reset_camera(L*0.1,L*0.9,L*0.1,L*0.9,-0.1,0.1)
    
    png_show(ren,filename='%s/init%04d'%(out_dir,run),width=800,height=800)	
    
    
    # <codecell>
    
    tau = monte_carlo_timestep(N_b,particles,U,params)
    f = open('%s/U%04d'%(out_dir,run), 'w')
    print tau
    f.write('%f\n'%tau)
    tau_new = monte_carlo_timestep(N_b,particles,U,params)
    while abs(tau_new-tau)/tau_new > tau_s:
        print tau_new
        f.write('%f\n'%tau_new)
        tau = tau_new
        tau_new = monte_carlo_timestep(N_b,particles,U,params)
        
    v = particles.get_grid()
    glyph.modified()
    png_show(ren,filename='%s/final%04d'%(out_dir,run),width=800,height=800)	
    w = tvtk.XMLUnstructuredGridWriter(input=v, file_name='%s/final%04d.vtu'%(out_dir,run))
    w.write()
    f.close()
    
    
    
