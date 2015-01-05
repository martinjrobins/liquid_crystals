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
import numpy as np
from random import sample
from utilities import import_columns



assert len(sys.argv)==2
out_dir = sys.argv[1]

L = 1
sigma_s = 0.5
k = 3




LdG_solutions = glob.glob('%s/LdG_solution_*.vtu'%out_dir)

averaged_files = glob.glob('%s/vtkAveraged*.vtu'%out_dir)
averaged_files.sort()

trajectory_files = glob.glob('%s/trajectories*.npy'%out_dir)

plt.figure(figsize=(6,4.5))
plt.xlabel('x')
plt.ylabel('y')
for filename in trajectory_files:
    print 'doing ',filename
    data = np.load(filename)
    #with np.load(filename) as data:
    plt.clf()
    if (data.shape[0] > 5):
        indicies = sample(range(data.shape[0]),5)
    else:
        indicies = range(data.shape[0])

    for i in indicies:
        plt.plot(data[i,0,:],data[i,1,:])
    plt.savefig('%s.pdf'%os.path.splitext(filename)[0])
        

if len(averaged_files)==0:
    averaged_files = glob.glob('%s/finalAveraged*.vtu'%out_dir)
    averaged_files.sort()
    #L = 1.0
    #sigma_s = L/200.0
    #k = 3
    
else:
    U_data = import_columns('%s/U'%out_dir)
    plt.figure(figsize=(6,4.5))
    plt.plot(U_data[0]+1,U_data[1], label="U")
    plt.xlabel("$N \times 10^5$")
    plt.ylabel("$U$")
    plt.savefig('%s/U_plot.pdf'%out_dir)
        

batch_files = glob.glob('%s/vtkBatch*.vtu'%out_dir) + glob.glob('%s/vtkInit.vtu'%out_dir)
batch_files.sort()
print batch_files
if len(batch_files)==0:
    batch_files = glob.glob('%s/finalBatch*.vtu'%out_dir)
    batch_files.sort()


#
# Common
#

lut = tvtk.LookupTable(hue_range=[0.66667,0.0],range=[0,1]);
lut.build()

lut2 = tvtk.LookupTable(hue_range=[0.66667,0.0],range=[0,0.1]);
lut2.build()

domain_source = tvtk.CubeSource(x_length=L,y_length=L,z_length=0,center=[L/2.0,L/2.0,0])
aDomainSource = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=domain_source.output),property=tvtk.Property(representation='w'))
aScalarBar = tvtk.ScalarBarActor(lookup_table=lut,title="s",position=[0.87,0.1],label_format='%-#2.1g',maximum_width_in_pixels=100)
aScalarBar2 = tvtk.ScalarBarActor(lookup_table=lut2,title="Q Var",position=[0.87,0.1],label_format='%-#2.1g',maximum_width_in_pixels=100)


transform = tvtk.Transform()
tvtk.to_vtk(transform).RotateZ(90)
cylinder_source = tvtk.TransformPolyDataFilter(input=tvtk.CylinderSource(height=sigma_s*2,radius=sigma_s/5).output,transform=transform)



#
# Averaged Visualisation
#
if len(averaged_files)>0:
    r = tvtk.XMLUnstructuredGridReader(file_name=averaged_files[0])
    calc1 = tvtk.ArrayCalculator(input=r.output,result_array_name="s",function="mag(averaged_orientation)",replacement_value=0, replace_invalid_values=True)
    calc1.add_vector_variable('averaged_orientation','averaged_orientation')
    calc2 = tvtk.ArrayCalculator(input=calc1.output,result_array_name="n1",function="sqrt(ao_X/(2*s) + 0.5)",replacement_value=0,replace_invalid_values=True)
    calc2.add_scalar_variable('ao_X','averaged_orientation',0)
    calc2.add_scalar_variable('s','s')
    calc3 = tvtk.ArrayCalculator(input=calc2.output,result_array_name="n2",function="ao_Y/(2*s*n1)",replacement_value=1,replace_invalid_values=True)
    calc3.add_scalar_variable('ao_Y','averaged_orientation',1)
    calc3.add_scalar_variable('s','s')
    calc3.add_scalar_variable('n1','n1')
    calc4 = tvtk.ArrayCalculator(input=calc3.output,result_array_name="n",function="n1*iHat + n2*jHat",replace_invalid_values=True)
    calc4.add_scalar_variable('n1','n1')
    calc4.add_scalar_variable('n2','n2')
    calc9 = tvtk.ArrayCalculator(input=r.output,result_array_name="Q11_var",function="variance_orientation_X/(n-1)",replacement_value=0, replace_invalid_values=True)
    calc9.add_scalar_variable('n','number_of_moves')
    calc9.add_scalar_variable('variance_orientation_X','variance_orientation',0)
    calc10 = tvtk.ArrayCalculator(input=calc9.output,result_array_name="Q12_var",function="variance_orientation_Y/(n-1)",replacement_value=0, replace_invalid_values=True)
    calc10.add_scalar_variable('n','number_of_moves')
    calc10.add_scalar_variable('variance_orientation_Y','variance_orientation',1)
    delaunay2 = tvtk.Delaunay2D(input=calc10.output)  
    delaunay2mapper = tvtk.PolyDataMapper(input=delaunay2.output,lookup_table=lut2)
    aDelaunay2 = tvtk.Actor(mapper=delaunay2mapper)
    ren4 = tvtk.Renderer()
    ren4.add_actor(aDelaunay2)
    ren4.add_actor(aScalarBar2)

    ren4.reset_camera(L*0.2,L,L*0.1,L*0.9,-0.1,0.1)
    
    renderWindow4 = tvtk.RenderWindow(off_screen_rendering=True,size=[850,800])
    renderWindow4.add_renderer(ren4)
    renderWindow4.render()
     
    windowToImageFilter4 = tvtk.WindowToImageFilter(input=renderWindow4)
    windowToImageFilter4.update()
    

    calc4.update()
    calc4.output.point_data.set_active_vectors('n')
    
    glyph = tvtk.Glyph3D(source=cylinder_source.output,input=calc4.output,scale_factor=1,vector_mode='use_vector',scaling=False,orient=True)
     
    calc1.update()
    calc1.output.point_data.set_active_scalars('s')
    delaunay = tvtk.Delaunay2D(input=calc1.output)
    
    aDelaunay = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=delaunay.output,lookup_table=lut))
    aGlyph = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=glyph.output, lookup_table=lut, scalar_visibility=False),property=tvtk.Property(color=(0.5,0.5,0.5)))
    
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

if len(batch_files) > 0:
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
    
    
#
# LdG Visualisation
#
if len(LdG_solutions) > 0:
    r3 = tvtk.XMLUnstructuredGridReader(file_name=LdG_solutions[0])
    
    calc5 = tvtk.ArrayCalculator(input=r3.output,result_array_name="s",function="mag(u)",replacement_value=0, replace_invalid_values=True)
    calc5.add_vector_variable('u','u')
    calc6 = tvtk.ArrayCalculator(input=calc5.output,result_array_name="n1",function="sqrt(ao_X/(2*s) + 0.5)",replacement_value=0,replace_invalid_values=True)
    calc6.add_scalar_variable('ao_X','u',0)
    calc6.add_scalar_variable('s','s')
    calc7 = tvtk.ArrayCalculator(input=calc6.output,result_array_name="n2",function="ao_Y/(2*s*n1)",replacement_value=1,replace_invalid_values=True)
    calc7.add_scalar_variable('ao_Y','u',1)
    calc7.add_scalar_variable('s','s')
    calc7.add_scalar_variable('n1','n1')
    calc8 = tvtk.ArrayCalculator(input=calc7.output,result_array_name="n",function="n1*iHat + n2*jHat",replace_invalid_values=True)
    calc8.add_scalar_variable('n1','n1')
    calc8.add_scalar_variable('n2','n2')
    
    calc8.update()
    calc8.output.point_data.set_active_vectors('n')
    glyph3 = tvtk.Glyph3D(source=cylinder_source.output,input=calc8.output,scale_factor=1,vector_mode='use_vector',scaling=False,orient=True)
    aGlyph3 = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=glyph3.output, lookup_table=lut, scalar_visibility=False),property=tvtk.Property(color=(0.5,0.5,0.5)))

    
    umesh = tvtk.GeometryFilter(input=calc5.output)
    umesh.update()
    umesh.output.point_data.set_active_scalars('s')
    aMesh = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=umesh.output,lookup_table=lut))

    ren3 = tvtk.Renderer()
    ren3.add_actor(aGlyph3)
    ren3.add_actor(aDomainSource)
    ren3.add_actor(aMesh)
    ren3.add_actor(aScalarBar)
    ren3.reset_camera(L*0.2,L,L*0.1,L*0.9,-0.1,0.1)
    
    renderWindow3 = tvtk.RenderWindow(off_screen_rendering=True,size=[850,800])
    renderWindow3.add_renderer(ren3)
    renderWindow3.render()
     
    windowToImageFilter3 = tvtk.WindowToImageFilter(input=renderWindow3)
    windowToImageFilter3.update()

     
for filename in averaged_files:
    print 'doing file ',filename
    r.file_name =filename
    r.modified()
    ren.render()
    windowToImageFilter.modified()
    windowToImageFilter.update()
    writer = tvtk.PNGWriter(file_name='%s.png'%os.path.splitext(filename)[0],input=windowToImageFilter.output)
    writer.write()
    
    calc10.update()
    calc10.output.point_data.set_active_scalars('Q11_var')
    delaunay2mapper.scalar_range = calc10.output.point_data.get_array('Q11_var').range
    delaunay2mapper.modified()
    ren4.render()
    windowToImageFilter4.modified()
    windowToImageFilter4.update()
    writer = tvtk.PNGWriter(file_name='%s_Q11_var.png'%os.path.splitext(filename)[0],input=windowToImageFilter4.output)
    writer.write()
    
    calc10.output.point_data.set_active_scalars('Q12_var')
    delaunay2mapper.scalar_range = calc10.output.point_data.get_array('Q12_var').range
    delaunay2mapper.modified()
    ren4.render()
    windowToImageFilter4.modified()
    windowToImageFilter4.update()
    writer = tvtk.PNGWriter(file_name='%s_Q12_var.png'%os.path.splitext(filename)[0],input=windowToImageFilter4.output)
    writer.write()
    
for filename in LdG_solutions:
    print 'doing file ',filename
    r3.file_name = filename
    r3.modified()
    ren3.render()
    windowToImageFilter3.modified()
    windowToImageFilter3.update()
    writer = tvtk.PNGWriter(file_name='%s.png'%os.path.splitext(filename)[0],input=windowToImageFilter3.output)
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
