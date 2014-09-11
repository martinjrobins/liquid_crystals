from particleSimulation import *
from tvtk.api import tvtk
from random import uniform
from math import sqrt

L = 1.0
k = 3.0
sigma_s = L/50.0

params = Params()
params['Dtrans'] = 1.0
params['Drot'] = 10000.0
params['Temp'] = 0.005
params['dt'] = 0.0000001
params['L'] = L

particles = Particles()
for i in range(100):
    p = Particle()
    p.position = Vect3d(uniform(0,L),uniform(0,L),0)
    u = [uniform(0,L),uniform(0,L),0]
    u = [i/sqrt(u[0]**2+u[1]**2) for i in u]
    p.orientation = Vect3d(u)
    particles.append(p)
 

U = GayBernePotential(sigma_s=sigma_s,k=k,kdash=1.0/5.0,mu=2,nu=1,epsilon_0=0.01)

monte_carlo_timestep(10,particles,U,params)

v = particles.get_grid()
v._get_point_data().set_active_normals('orientation')
arrow = tvtk.SphereSource()
glyph = tvtk.Glyph3D(source=arrow.output,input=v,scale_factor=sigma_s*k,vector_mode=True,orient=True)
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
