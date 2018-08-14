import vtk
from vtk.util import numpy_support as VN
import numpy as np
#from numba import jit

def vtiread(filename):
    imageData = vtk.vtkImageData()
    reader    = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()

    data    = reader.GetOutput()
    spacing = np.array(data.GetSpacing())
    extent  = np.array(data.GetExtent())
    dim     = data.GetDimensions()
    dim     = np.array([i-1 for i in dim])

    u = VN.vtk_to_numpy(data.GetCellData().GetArray('facedata:01:01'))
    v = VN.vtk_to_numpy(data.GetCellData().GetArray('facedata:01:02'))
    w = VN.vtk_to_numpy(data.GetCellData().GetArray('facedata:01:03'))

    #u = u.reshape(vec,order='F')
    #v = v.reshape(vec,order='F')
    #w = w.reshape(vec,order='F')

    #x = np.zeros(data.GetNumberOfPoints())
    #y = np.zeros(data.GetNumberOfPoints())
    #z = np.zeros(data.GetNumberOfPoints())


    origin = np.array(data.GetPoint(0))

    #z,y,x = np.meshgrid(np.arange(dim[0]), np.arange(dim[1]),
            #np.arange(dim[2]))

    #x = x * spacing[0] + origin[0]
    #y = y * spacing[1] + origin[1]
    #z = z * spacing[2] + origin[2]

    vel = np.stack((u,v,w), axis = 1)
    #coord = np.stack((x.reshape(-1),y.reshape(-1),z.reshape(-1)), axis = 1)

    #x = x.reshape(dim,order='F')
    #y = y.reshape(dim,order='F')
    #z = z.reshape(dim,order='F')

    return vel, origin, spacing, dim
