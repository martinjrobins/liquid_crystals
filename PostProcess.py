# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 00:48:17 2014

@author: mrobins
"""

import glob
from tvtk.api import tvtk
from png_show import png_show


out_dir = '/home/mrobins/tmp/HGO/new_averaging'
sigma_s = 0.5
k = 3
L = 25

batch_files = glob.glob('%s/vtkAveraged*.vtu'%out_dir)
batch_files.sort()
batch_files = batch_files[1:-1]
print batch_files
for filename in batch_files:
    print 'doing file ',filename
    r = tvtk.XMLUnstructuredGridReader(file_name=filename)
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
    calc4.update()
    calc4.output.point_data.set_active_scalars('s')
    calc4.output.point_data.set_active_vectors('n')

    
    ellipse = tvtk.ParametricEllipsoid(x_radius=sigma_s*k/2.0,y_radius=sigma_s/2.0,z_radius=sigma_s/2.0)
    ellipse_source = tvtk.ParametricFunctionSource(parametric_function=ellipse,u_resolution=8,v_resolution=8,w_resolution=8)
    domain_source = tvtk.CubeSource(x_length=L,y_length=L,z_length=0,center=[L/2.0,L/2.0,0])
    m2 = tvtk.PolyDataMapper(input=domain_source.output)
    p = tvtk.Property(representation='w')
    a2 = tvtk.Actor(mapper=m2,property=p)

    glyph2 = tvtk.Glyph3D(source=ellipse_source.output,input=calc4.output,scale_factor=1,vector_mode=False,scaling=False,orient=True)
    #glyph2.set_input_array_to_process(1,0,0,0,'n')	# vectors
   # glyph2.set_input_array_to_process(3,0,0,0,'s')	# colors
    glyph2.update()
    
    delaunay = tvtk.Delaunay2D(input=calc1.output)
    delaunay.update()
    delaunay.output.point_data.set_active_scalars('s')


    lut = tvtk.LookupTable(hue_range=[0.66667,0.0],range=[0,1]);
    lut.build()

    mDelaunay = tvtk.PolyDataMapper(input=delaunay.output,lookup_table=lut)    
    m3 = tvtk.PolyDataMapper(input=glyph2.output, lookup_table=lut)
    a4 = tvtk.ScalarBarActor(lookup_table=lut,title="s")
    
    a3 = tvtk.Actor(mapper=m3)
    ren2 = tvtk.Renderer()
    ren2.add_actor(a3)
    ren2.add_actor(a2)
    ren2.add_actor(a4)
    ren2.add_actor(tvtk.Actor(mapper=mDelaunay))
    ren2.reset_camera(L*0.1,L*0.9,L*0.1,L*0.9,-0.1,0.1)
    renderWindow = tvtk.RenderWindow(size=[800,600])
    renderWindow.add_renderer(ren2)
    renderWindowInteractor = tvtk.RenderWindowInteractor(render_window = renderWindow)
    renderWindowInteractor.start()
    break