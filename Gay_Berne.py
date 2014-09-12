from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt,pi,cos,sin

L = 50.0
k = 3.0
sigma_s = 1.0
rho = 0.7
rot_step = 2*pi/300
diff_step = sigma_s/50
dt = 0.004
T = 3.2
#vol = (4.0/3.0)*pi*k*sigma_s**3
area = (1.0/4.0)*pi*k*sigma_s**2
N = int(rho*L**2/area)

params = Params()
params['Dtrans'] = diff_step**2/(2.0*dt)*100
params['Drot'] = rot_step**2/(2.0*dt)*100
params['Temp'] = T
params['dt'] = dt
params['L'] = L

particles = Particles()
for i in range(N):
    p = Particle()
    p.position = Vect3d(uniform(0,L),uniform(0,L),0)
    theta = uniform(0,2*pi)
    u = [cos(theta),sin(theta),0]
    p.orientation = Vect3d(u)
    particles.append(p)
 
U_hgo = HGOPotential(sigma_s=sigma_s,k=k)
U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=1,nu=3,epsilon_0=1)

monte_carlo_timestep(1000,particles,U_hgo,params)

v = particles.get_grid()
v._get_point_data().set_active_normals('orientation')
ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
glyph = tvtk.Glyph3D(source=ellipse_source.output,input=v,scale_factor=1,vector_mode=True,orient=True)

m = tvtk.PolyDataMapper(input=glyph.output)
a = tvtk.Actor(mapper=m)
ren = tvtk.Renderer()
ren.add_actor(a)
renWin = tvtk.RenderWindow(size=[800,600])
renWin.add_renderer(ren)
renWin.render()

params['Dtrans'] = diff_step**2/(2.0*dt)
params['Drot'] = rot_step**2/(2.0*dt)
for i in range(1000):
    print i
    monte_carlo_timestep(100,particles,U,params)
    v = particles.get_grid()
    glyph.modified()
    renWin.render()
    

#renWinInt = tvtk.RenderWindowInteractor(render_window=renWin)
#renWinInt.start()
