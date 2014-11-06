# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 17:14:17 2014

@author: robinsonm
"""

import sys
from tvtk.api import tvtk
from particleSimulation import *

assert len(sys.argv)==4

out_dir = sys.argv[1]
fromi = int(sys.argv[2])
toi = int(sys.argv[3])

print "out_dir = ",out_dir," averaging from i = ",fromi," to ",toi

particles_sum = Particles()
particles = Particles()

for i in range(fromi,toi):
    filename = '%s/vtkAveraged%04d.vtu'%(out_dir,i)
    print 'reading from file ',filename

    r = tvtk.XMLUnstructuredGridReader(file_name=filename)
    r.update()
    grid = tvtk.to_vtk(r.output)
    particles.copy_from_vtk_grid(grid)
    
    if i == fromi:
        particles_sum.copy_from_vtk_grid(tvtk.to_vtk(r.output))
    else:
        for p_sum,p in zip(particles_sum,particles):
            p_sum.averaged_orientation = p_sum.averaged_orientation + p.averaged_orientation

for p_sum in particles_sum:
    p_sum.averaged_orientation = p_sum.averaged_orientation / (toi-fromi)
        
w = tvtk.XMLUnstructuredGridWriter(input=particles_sum.get_grid(), file_name='%s/vtkAveragedAverage%04d_%04d.vtu'%(out_dir,fromi,toi))
w.write()