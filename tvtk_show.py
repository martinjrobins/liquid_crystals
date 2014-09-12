# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 12:31:25 2014

@author: robinsonm
"""

from IPython.display import Image
from tvtk.api import tvtk
def tvtk_show(renderer, width=400, height=300):
    """
    Takes vtkRenderer instance and returns an IPython Image with the rendering.
    """
    renderWindow = tvtk.RenderWindow(off_screen_rendering=True,size=[width,height])
    renderWindow.add_renderer(renderer)
    renderWindow.render()
     
    windowToImageFilter = tvtk.WindowToImageFilter(input=renderWindow)
    windowToImageFilter.update()
     
    writer = tvtk.PNGWriter(write_to_memory=True,input=windowToImageFilter.output)
    writer.write()
    data = str(buffer(tvtk.to_vtk(writer).GetResult()))
    return Image(data)