from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi

L = 50.0
k = 3.0
sigma_s = 1.0
rho = 0.5
rot_step = 2*pi/10
diff_step = sigma_s/5
dt = 0.004
T = 3.2
#vol = (4.0/3.0)*pi*k*sigma_s**3
area = (1.0/4.0)*pi*k*sigma_s**2
N = int(rho*L**2/area)

params = Params()
params['Dtrans'] = diff_step**2/(2.0*dt)
params['Drot'] = rot_step**2/(2.0*dt)
params['Temp'] = T
params['dt'] = dt
params['L'] = L

particles = Particles()
for i in range(N):
    p = Particle()
    p.position = Vect3d(uniform(0,L),uniform(0,L),0)
    u = [uniform(0,L),uniform(0,L),0]
    u = [i/sqrt(u[0]**2+u[1]**2) for i in u]
    p.orientation = Vect3d(u)
    particles.append(p)
 

U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)

monte_carlo_timestep(10,particles,U,params)

v = particles.get_grid()
v._get_point_data().set_active_normals('orientation')
#arrow = tvtk.ArrowSource(tip_length=sigma_s*k)
sphere = tvtk.SphereSource(radius=sigma_s)
ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
glyph = tvtk.Glyph3D(source=ellipse_source.output,input=v,scale_factor=1,vector_mode=True,orient=True)
#glyph2 = tvtk.Glyph3D(source=ellipse_source.output,input=v,scale_factor=1,vector_mode=True,orient=True)

m = tvtk.PolyDataMapper(input=glyph.output)
a = tvtk.Actor(mapper=m)
ren = tvtk.Renderer()
ren.add_actor(a)
renWin = tvtk.RenderWindow(size=[800,600])
renWin.add_renderer(ren)
renWin.render()

for i in range(1000):
    print i
    monte_carlo_timestep(10,particles,U,params)
    v = particles.get_grid()
    #glyph.input=v
    glyph.modified()
    renWin.render()
    

#renWinInt = tvtk.RenderWindowInteractor(render_window=renWin)
#renWinInt.start()
