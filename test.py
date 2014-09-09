from particleSimulation import *
from tvtk.api import tvtk

L = 1.0

sim = ParticleSimulation(Drot=1.0,Dtrans=1.0,Temp=1.0,dt=0.001,L=L)

U = GayBernePotential(sigma_s=L/100.0,k=3,kdash=1.0/5.0,mu=2,nu=1,epsilon_0=1.0)

sim.add_particles(100)
sim.monte_carlo_timestep(1,U)
v = sim.particles.get_grid()
