# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 00:48:17 2014

@author: mrobins
"""

import glob
from tvtk.api import tvtk
import os
import sys
import matplotlib.pyplot as plt
from utilities import import_columns



assert len(sys.argv)==2
out_dir = sys.argv[1]

sigma_s = 0.5
k = 3
L = 25


U_data = import_columns('%s/U0000'%out_dir)
plt.figure(figsize=(6,4.5))
plt.plot(U_data[0]+1,U_data[1], label="U")
plt.xlabel("$N \times 10^5$")
plt.ylabel("$U$")
plt.savefig('%s/U_plot.pdf'%out_dir)


averaged_files = glob.glob('%s/vtkAveraged*.vtu'%out_dir)
averaged_files.sort()

batch_files = glob.glob('%s/vtkBatch*.vtu'%out_dir)
batch_files.sort()


#
# Averaged Visualisation
#

r = tvtk.XMLUnstructuredGridReader(file_name=averaged_files[0])
calc1 = tvtk.ArrayCalculator(input=r.output,result_array_name="s",function="mag(averaged_orientation)",replace_invalid_values=True)
calc1.add_vector_variable('averaged_orientation','averaged_orientation')
calc2 = tvtk.ArrayCalculator(input=calc1.output,result_array_name="n1",function="sqrt(ao_X/(2*s) + 0.5)",replace_invalid_values=True)
calc2.add_scalar_variable('ao_X','averaged_orientation',0)
calc2.add_scalar_variable('s','s')
calc3 = tvtk.ArrayCalculator(input=calc2.output,result_array_name="n2",function="ao_Y/(2*s*n1)",replace_invalid_values=True)
calc3.add_scalar_variable('ao_Y','averaged_orientation',1)
calc3.add_scalar_variable('s','s')
calc3.add_scalar_variable('n1','n1')
calc4 = tvtk.ArrayCalculator(input=calc3.output,result_array_name="n",function="n1*iHat + n2*jHat",replace_invalid_values=True)
calc4.add_scalar_variable('n1','n1')
calc4.add_scalar_variable('n2','n2')
   

#ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
#ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
domain_source = tvtk.CubeSource(x_length=L,y_length=L,z_length=0,center=[L/2.0,L/2.0,0])

calc4.update()
calc4.output.point_data.set_active_vectors('n')

transform = tvtk.Transform()
tvtk.to_vtk(transform).RotateZ(90)
cylinder_source = tvtk.TransformPolyDataFilter(input=tvtk.CylinderSource(height=1,radius=0.1).output,transform=transform)
glyph = tvtk.Glyph3D(source=cylinder_source.output,input=calc4.output,scale_factor=1,vector_mode='use_vector',scaling=False,orient=True)
 
calc1.update()
calc1.output.point_data.set_active_scalars('s')
delaunay = tvtk.Delaunay2D(input=calc1.output)

lut = tvtk.LookupTable(hue_range=[0.66667,0.0],range=[0,1]);
lut.build()

aDelaunay = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=delaunay.output,lookup_table=lut))
aScalarBar = tvtk.ScalarBarActor(lookup_table=lut,title="s",position=[0.87,0.1],label_format='%-#2.1g',maximum_width_in_pixels=100)
aGlyph = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=glyph.output, lookup_table=lut, scalar_visibility=False),property=tvtk.Property(color=(0.5,0.5,0.5)))
aDomainSource = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=domain_source.output),property=tvtk.Property(representation='w'))

ren = tvtk.Renderer()
ren.add_actor(aDelaunay)
ren.add_actor(aScalarBar)
ren.add_actor(aGlyph)
#ren.add_actor(aDomainSource)

ren.reset_camera(L*0.2,L,L*0.1,L*0.9,-0.1,0.1)

renderWindow = tvtk.RenderWindow(off_screen_rendering=True,size=[850,800])
renderWindow.add_renderer(ren)
renderWindow.render()
 
windowToImageFilter = tvtk.WindowToImageFilter(input=renderWindow)
windowToImageFilter.update()


#
# Batch Visualisation
#

r2 = tvtk.XMLUnstructuredGridReader(file_name=batch_files[0])

ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)

r2.update()
r2.output.point_data.set_active_vectors('orientation')
glyph2 = tvtk.Glyph3D(source=ellipse_source.output,input=r2.output,scale_factor=1,vector_mode='use_vector',scaling=False,orient=True)

aGlyph2 = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=glyph2.output, lookup_table=lut, scalar_visibility=False),property=tvtk.Property(color=(1.0,1.0,1.0)))

ren2 = tvtk.Renderer()
ren2.add_actor(aGlyph2)
ren2.add_actor(aDomainSource)
ren2.reset_camera(L*0.1,L*0.9,L*0.1,L*0.9,-0.1,0.1)

renderWindow2 = tvtk.RenderWindow(off_screen_rendering=True,size=[800,800])
renderWindow2.add_renderer(ren2)
renderWindow2.render()
 
windowToImageFilter2 = tvtk.WindowToImageFilter(input=renderWindow2)
windowToImageFilter2.update()

     
for filename in averaged_files:
    print 'doing file ',filename
    r.file_name =filename
    r.modified()
    ren.render()
    windowToImageFilter.modified()
    windowToImageFilter.update()
    writer = tvtk.PNGWriter(file_name='%s.png'%os.path.splitext(filename)[0],input=windowToImageFilter.output)
    writer.write()
    
for filename in batch_files:
    print 'doing file ',filename
    r2.file_name =filename
    r2.modified()
    r2.update()
    r2.output.point_data.set_active_vectors('orientation')
    ren2.render()
    windowToImageFilter2.modified()
    windowToImageFilter2.update()
    writer = tvtk.PNGWriter(file_name='%s.png'%os.path.splitext(filename)[0],input=windowToImageFilter2.output)
    writer.write()
