import os
import numpy as np

from read_node import *

def appendSpherical(coord):
    ptsnew = np.hstack( (coord, np.zeros(coord.shape)) )
    xy = coord[:,0]**2 + coord[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + coord[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), coord[:,2]) # for elevation angle defined from Z-axis down
    ptsnew[:,5] = np.arctan2(coord[:,1], coord[:,0])
    return ptsnew

center = np.array([18.8495559215, 18.8495559215, 18.8495559215])

# get a sample
filename = "flow/node-0000000001.vti"
vel, origin, spacing, dim = vtiread(filename)

print(vel)
coord = coord - center
coord_6 = appendSpherical(coord)

print(coord_6[:,3])
print(coord_6.shape)


# Stream

directory = os.getcwd()
directory += "/flow"

### Counting how many files
f_count = len([f for f in os.listdir(directory)
     if f.endswith('.vti') and os.path.isfile(os.path.join(directory, f))])

print(f_count)


coord_all = np.zeros(coord_6.shape)
vel_all   = np.zeros(vel.shape)
print(coord_6.shape)

### Everything together

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename.endswith(".vti"):
        vti_name = os.path.join(directory, filename)

        coord, vel = vtiread(vti_name)
        coord      = coord - center
        coord_6    = appendSpherical(coord)

        coord_all  = np.hstack((coord_all, coord_6))
        vel_all    = np.hstack((vel_all, vel))
        continue

    else:
        continue



