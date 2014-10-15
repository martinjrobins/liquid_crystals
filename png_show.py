# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:18:23 2014

@author: robinsonm
"""

from IPython.display import Image
from tvtk.api import tvtk
def png_show(renderer, filename='test.png', width=400, height=300):
    """
    Takes vtkRenderer instance and outputs a png file
    """
    renderWindow = tvtk.RenderWindow(off_screen_rendering=True,size=[width,height])
    renderWindow.add_renderer(renderer)
    renderWindow.render()
     
    windowToImageFilter = tvtk.WindowToImageFilter(input=renderWindow)
    windowToImageFilter.update()
     
    writer = tvtk.PNGWriter(file_name=filename,input=windowToImageFilter.output)
    writer.write()