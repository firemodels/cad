# EVM, March 2022
# This script creates a binary bulk density file that can be read by FDS
# It is an example case, but can be repurposed to handle 
# specific tree/canopy which the user might have

import numpy as np
from scipy.io import FortranFile

# Create a binary file to hold data
f = FortranFile('canopy_foliage.bdf', 'w')

###############################################################################
# This section loads a sample 3D array of bulk density data, 350 kg total
# In reality the user would provide this data from another source

bd_data=np.load('./bd_data.npy')

# array/voxel resolution (m)
dx=dy=dz=1
# size of array
nx,ny,nz=bd_data.shape

# meshgrid of x,y,z voxel centers
xv=np.linspace(dx/2,dx*(nx-1/2),nx)
yv=np.linspace(dy/2,dy*(ny-1/2),ny)
zv=np.linspace(dz/2,dz*(nz-1/2),nz)
xx,yy,zz = np.meshgrid(xv,yv,zv,indexing='ij')

###############################################################################

# relevant voxel centers
xv=xx[bd_data>0]
yv=yy[bd_data>0]
zv=zz[bd_data>0]

# write out global bounding voxel faces
VXMIN=min(xv)-dx/2
VXMAX=max(xv)+dx/2
VYMIN=min(yv)-dy/2
VYMAX=max(yv)+dy/2
VZMIN=min(zv)-dz/2
VZMAX=max(zv)+dz/2
f.write_record(np.array([VXMIN,VXMAX,VYMIN,VYMAX,VZMIN,VZMAX],dtype=np.float64))

# write out voxel resolution
f.write_record(np.array([dx,dy,dz],dtype=np.float64))

# write out number of relevant (non-zero) voxels
bd_data=bd_data[bd_data>0]
nvox=bd_data.shape[0]
f.write_record(np.array(nvox,dtype=np.int32))

# write out center and bulk density for each non-zero voxel
for (vxc,vyc,vzc,bd) in zip(xv,yv,zv,bd_data):
    xyz=[vxc,vyc,vzc]
    f.write_record(np.array(xyz,dtype=np.float64))
    f.write_record(np.array(bd,dtype=np.float64))
    
f.close()


    

    

    

    

    


