# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 17:09:56 2014

@author: robinsonm
"""

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import array
import vtk



def import_columns(filename):
    reader = csv.reader(open(filename), delimiter=' ',skipinitialspace=True)
    first_line = reader.next()
    print first_line
    decrement = 0
    if (first_line[0][0]=='#'):
        decrement = decrement+1
    if (first_line[-1]==''):
        decrement = decrement+1
        
    n = len(first_line)-decrement
    columns = [array.array('f') for _ in xrange(n)]
    
    if first_line[0][0]!='#':
        for i in xrange(n):
                columns[i].append(float(first_line[i]))
    for line in reader:
        if line[0][0]=='#':
            continue
        for i in xrange(n):
            columns[i].append(float(line[i]))
    for i in xrange(n):
        columns[i] = np.array(columns[i])                                      
    return columns

def mesh(x,y,z):
    nx = 100
    ny = 100
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    xi, yi = np.meshgrid(xi, yi) 
    
    x_new = (x - xmin) / (xmax - xmin)
    xi_new = (xi - xmin) / (xmax - xmin)
    y_new = (y - ymin) / (ymax - ymin)
    yi_new = (yi - ymin) / (ymax - ymin)

    print x_new.shape
    print y_new.shape
    print z.shape
    print xi_new.shape
    print yi_new.shape
    zi = plt.mlab.griddata(x_new, y_new, z, xi_new, yi_new)

    return xi,yi,zi

    
def volume_plot(input_data_matrix):
    nx,ny,nz = input_data_matrix.shape
    min_data = input_data_matrix.min()
    max_data = input_data_matrix.max()

    print 'volume plot of data matrix with shape ',input_data_matrix.shape,' and (min,max) = (',min_data,',',max_data,")"

    data_matrix = np.uint8(255.0*(input_data_matrix-min_data)/(max_data-min_data))
#    data_matrix[0:35, 0:35, 0:35] = 50
#    data_matrix[25:55, 25:55, 25:55] = 100
#    data_matrix[45:69, 45:69, 45:69] = 150
    
    nx,ny,nz = data_matrix.shape
    min_data = data_matrix.min()
    max_data = data_matrix.max()
    
    #print 'volume plot of data matrix with shape ',data_matrix.shape,' and (min,max) = (',min_data,',',max_data,")"
    # For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
    # imports raw data and stores it. 
    dataImporter = vtk.vtkImageImport()
    # The preaviusly created array is converted to a string of chars and imported.
    data_string = data_matrix.tostring()
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToUnsignedChar()
    # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
    # must be told this is the case.
    dataImporter.SetNumberOfScalarComponents(1)
    # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
    # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
    # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
    # VTK complains if not both are used.
    dataImporter.SetDataExtent(0, nx-1, 0, ny-1, 0, nz-1)
    dataImporter.SetWholeExtent(0, nx-1, 0, ny-1, 0, nz-1)
    
    # The following class is used to store transparencyv-values for later retrival. In our case, we want the value 0 to be
    # completly opaque whereas the three different cubes are given different transperancy-values to show how it works.
    alphaChannelFunc = vtk.vtkPiecewiseFunction()
    alphaChannelFunc.AddPoint(min_data, 0.0)
    alphaChannelFunc.AddPoint(min_data+1, 0.2)
    #alphaChannelFunc.AddPoint(0.5*(max_data+min_data), 0.1)
    #alphaChannelFunc.AddPoint(max_data, 1.0)
    
#    alphaChannelFunc.AddPoint(0, 0.0)
#    alphaChannelFunc.AddPoint(50, 0.05)
#    alphaChannelFunc.AddPoint(100, 0.1)
#    alphaChannelFunc.AddPoint(150, 0.2)

    
    # This class stores color data and can create color tables from a few color points. For this demo, we want the three cubes
    # to be of the colors red green and blue.
    colorFunc = vtk.vtkColorTransferFunction()
    print min_data,max_data
    colorFunc.AddRGBPoint(min_data, 1.0, 0.0, 0.0)
    colorFunc.AddRGBPoint(max_data, 0.0, 0.0, 1.0)
    
    # The preavius two classes stored properties. Because we want to apply these properties to the volume we want to render,
    # we have to store them in a class that stores volume prpoperties.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(colorFunc)
    volumeProperty.SetScalarOpacity(alphaChannelFunc)
    
    
    # This class describes how the volume is rendered (through ray tracing).
    compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
    # We can finally create our volume. We also have to specify the data for it, as well as how the data will be rendered.
    volumeMapper = vtk.vtkVolumeRayCastMapper()
    volumeMapper.SetVolumeRayCastFunction(compositeFunction)
    volumeMapper.SetInputConnection(dataImporter.GetOutputPort())
    
    # The class vtkVolume is used to pair the preaviusly declared volume as well as the properties to be used when rendering that volume.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)
    
    # With almost everything else ready, its time to initialize the renderer and window, as well as creating a method for exiting the application
    renderer = vtk.vtkRenderer()
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)
    
    # We add the volume to the renderer ...
    renderer.AddVolume(volume)
    # ... set background color to white ...
    renderer.SetBackground(1, 1, 1)
    # ... and set window size.
    renderWin.SetSize(400, 400)
    
    # A simple function to be called when the user decides to quit the application.
    def exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)
    
    # Tell the application to use the function as an exit check.
    renderWin.AddObserver("AbortCheckEvent", exitCheck)
    
    renderInteractor.Initialize()
    # Because nothing will be rendered without any input, we order the first render manually before control is handed over to the main-loop.
    renderWin.Render()
    renderInteractor.Start()

