import sys
import gridData 
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
import numpy as np

g=gridData.Grid('src/test/resources/mol1_HDON.dx')

T=[[-0.9730363128062556,-0.1469271768646341,-0.17779971501382005,-0.40339353445144477],
 [0.23053096053663238,-0.5945393991281284,-0.7703105731576206,3.039172174760854],
 [0.007470622093558343,-0.7905284989062399,0.6123796879572205,-5.433572743242558],
 [0.0,0.0,0.0,1.0] ]
T=np.array(T)
R=T[:3,:3]

# Sample data from DX file (replace this with your actual data)
dimensions = g.grid.shape
dx, dy, dz = dimensions
origin = g.origin
cell_spacing = g.delta

# Create original grid coordinates
x = np.arange(origin[0], origin[0] + dimensions[0] * cell_spacing[0], cell_spacing[0])
y = np.arange(origin[1], origin[1] + dimensions[1] * cell_spacing[1], cell_spacing[1])
z = np.arange(origin[2], origin[2] + dimensions[2] * cell_spacing[2], cell_spacing[2])

# Generate original grid in 3D space
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

# Convert the 3D grid to a 4D grid by adding a dimension with 1's
original_coordinates = np.stack((xx, yy, zz, np.ones_like(xx)), axis=-1)

# Flatten the original coordinates for matrix multiplication
flat_original_coords = original_coordinates.reshape((-1, 4)).T

# Apply the transformation matrix to the coordinates
transformed_coordinates = np.matmul(T.T, flat_original_coords).T

# Flatten the original density values for interpolation
flat_grid_values = g.grid.flatten()

transformed_density_values = griddata(
    (xx.flatten(), yy.flatten(), zz.flatten()),
    flat_grid_values,
    (transformed_coordinates[..., 0], transformed_coordinates[..., 1], transformed_coordinates[..., 2]),
    method='linear',
    fill_value=0.0  # Replace this with an appropriate fill value
)
new_vol = transformed_density_values.reshape(dimensions)
g.grid = new_vol
g.export('mol1_HDON_rotated.dx')


coords=np.meshgrid(np.arange(dimensions[0]),np.arange(dimensions[1]), np.arange(dimensions[2]))
xyz = np.vstack([coords[0].flatten(), coords[1].flatten(), coords[2].flatten(),np.ones(dx*dy*dz)])

# stack the meshgrid to position vectors, center them around 0 by substracting dim/2
xyz=np.vstack([coords[0].reshape(-1)-float(dx)/2,     # x coordinate, centered
                coords[1].reshape(-1)-float(dy)/2,     # y coordinate, centered
                coords[2].reshape(-1)-float(dz)/2] )    # z coordinate, centered
                # )    # 1 for homogeneous coordinates
               

transformed_xyz=np.matmul(T,xyz)
# new_xyz = transformed_xyz[:3,:].T
# extract coordinates, don't use transformed_xyz[3,:] that's the homogeneous coordinate, always 1
x=transformed_xyz[0,:]+float(dx)/2
y=transformed_xyz[1,:]+float(dy)/2
z=transformed_xyz[2,:]+float(dz)/2

# x=x.reshape(dimensions)
# y=y.reshape(dimensions)
# z=z.reshape(dimensions)
x=transformed_xyz[0,:].reshape(dimensions)
y=transformed_xyz[1,:].reshape(dimensions)
z=transformed_xyz[2,:].reshape(dimensions)

# the coordinate system seems to be strange, it has to be ordered like this
new_xyz=[x,y,z]

# sample
g=gridData.Grid('src/test/resources/mol1_HDON.dx')
#new_vol=scipy.ndimage.map_coordinates(g.grid,new_xyz)
transformed_density_values = griddata(
    (coords[0].flatten(), coords[1].flatten(), coords[2].flatten()),
    g.grid.flatten(),
    (transformed_xyz[0,:], transformed_xyz[1,:], transformed_xyz[2,:]),
    method='linear',
    fill_value=0.0  # Replace this with an appropriate fill value
)
new_vol = transformed_density_values.reshape(dimensions)
g.grid = new_vol
g.export('mol1_HDON_rotated.dx')