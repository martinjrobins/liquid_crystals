from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi,cos,sin

L = 50.0
rot_step = 10e-2
diff_step = 0
dt = 0.004
T = 0.05
N = 50

params = Params()
params['Dtrans'] = diff_step**2/(2.0*dt)
params['Drot'] = rot_step**2/(2.0*dt)
params['Temp'] = T
params['dt'] = dt
params['L'] = L

particles = Particles()
for i in range(N-1):
    for j in range(N-1): 
        p = Particle()
        p.position = Vect3d((i+1),(j+1),0)
        theta = uniform(0,2*pi)
        u = [cos(theta),sin(theta),0]
        p.orientation = Vect3d(u)
        particles.append(p)
 

U = LabwohlLasherPotential(epsilon=1,lattice_spacing=1)

monte_carlo_timestep(500,particles,U,params)

v = particles.get_grid()
v._get_point_data().set_active_normals('orientation')
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

renWin = tvtk.RenderWindow(size=[800,600])
renWin.add_renderer(ren)
renWin.render()

for i in range(100):
    print i
    monte_carlo_timestep(500,particles,U,params)
    v = particles.get_grid()
    glyph.modified()
    renWin.render()
    

#renWinInt = tvtk.RenderWindowInteractor(render_window=renWin)
#renWinInt.start()
