from osgeo import gdal
from gdalconst import *
import numpy as np
import struct
import scipy.ndimage.filters as filt

# Michael Christoffersen
# Python 2.7
# Uses GDAL, numpy, scipy
# Turns a DEM into an STL file with a base

dem_path = './megt_n_512.tif'
out_path = './megt_n_512.stl'
dsamp = 10 #factor for downsampling the DEM
vx = 50 #vertical exaggeration
volume = True #if false output is a surface, if true sides and a base are added
base = 4000 #this is elevation below the max height to put the base, in DEM units
smooth = True # Whether to smooth with a gaussian kernel
sigma = 3 # For smoothing filter

# Use this if the DEM has bad values stored for the pixel dimensions
# The value should be in the same units as the DEM heights 
force_dim = False
dim_val = 30 

### No Editing Below Here ###

# Cross product of two vectors
def xprd(v1, v2):
	i = v1[1]*v2[2] - v1[2]*v2[1]
	j = -(v1[0]*v2[2] - v1[2]*v2[0])
	k = v1[0]*v2[1] - v1[1]*v2[0]
	return [i,j,k]

# Read in the DEM
src_dem = gdal.Open(dem_path,GA_ReadOnly)
dem_data = np.array(src_dem.GetRasterBand(1).ReadAsArray()).astype('float32')

# Make sure elevations are > 0
if(dem_data.min() < 0):
	dem_data = dem_data + dem_data.min() + 1

# Subsample it and apply vertical exaggeration 
dem_data = dem_data[::dsamp, ::dsamp]*vx

if(smooth): # Filter the dem
	dem_data = filt.gaussian_filter(dem_data,sigma)

# Get DEM dimensions and pixel dimensions
y,x = dem_data.shape
print("{} x {}".format(x,y))
gt = src_dem.GetGeoTransform()

# This first part is for forcing pixel dimensions on poorly formed DEMs
# Otherwise, just read the geotransform
if(force_dim):
	xres = dim_val * dsamp
	yres = dim_val * dsamp
else:
	xres = gt[1]*dsamp
	yres = -gt[5]*dsamp
print("xres: {}".format(xres))
print("yres: {}".format(yres))


# Make an array of vectors along the x axis and y axis, pointing from point to adjacent point
yxvec = np.zeros((y-1,x))
yyvec = np.ones((y-1,x))*yres
yzvec = np.zeros((y-1,x))

xxvec = np.ones((y,x-1))*xres
xyvec = np.zeros((y,x-1))
xzvec = np.zeros((y,x-1))

for i in range(y-1):
	np.subtract(dem_data[i+1,:],dem_data[i,:],yzvec[i,:])

for i in range(x-1):
	np.subtract(dem_data[:,i+1],dem_data[:,i],xzvec[:,i])

dy = np.stack((yxvec,yyvec,yzvec),0)
dx = np.stack((xxvec,xyvec,xzvec),0)

# Start making STL file
outfl = open(out_path,'wb')

head = struct.pack('<'+'d'*10,0,0,0,0,0,0,0,0,0,0) # 80 byte padding at beginning
outfl.write(head)

# Number of facets for making a volume vs a surface
if(volume):
	ftot = 2*(x-1)*(y-1) + 4*(x-1) + 4*(y-1) + 2 # Total number of facets
else:
	ftot = (x-1)*(y-1)*2

print("{} facets".format(ftot))
leng = struct.pack('<I',ftot)
outfl.write(leng)

# Make the surface
# Making triangles and calculating normal vectors
for i in range(x-1):
	for j in range(y-1):
		ptl = dem_data[j,i]
		pbl = dem_data[j+1,i]
		ptr = dem_data[j,i+1]
		pbr = dem_data[j+1,i+1]
		vt = dx[:,j,i]
		vb = dx[:,j+1,i]
		vl = dy[:,j,i]
		vr = dy[:,j,i+1]
		xptl = xprd(vt,vl)
		xpbr = xprd(vb,vr)
		ftl = struct.pack('<' + 'f'*12 + 'H',
				xptl[0],xptl[1],xptl[2],
				i*xres, (j+1)*yres, pbl,
				i*xres, j*yres, ptl,
				(i+1)*xres, j*yres, ptr, 0)
		fbr = struct.pack('<' + 'f'*12 + 'H',
				xpbr[0],xpbr[1],xpbr[2],
				(i+1)*xres, (j+1)*yres, pbr,
				i*xres, (j+1)*yres, pbl,
				(i+1)*xres, j*yres, ptr, 0)
		#print(ptl,pbl,ptr,pbr)
		outfl.write(ftl)
		outfl.write(fbr)

# Add sides and a base if making a volume
if(volume):
	bz = dem_data.max() - base*vx

	# Add two sides
	for i in range(x-1):
		ptl = dem_data[0,i]
		ptr = dem_data[0,i+1]
		pbl = dem_data[y-1,i]
		pbr = dem_data[y-1,i+1]
		t1 = struct.pack('<' + 'f'*12 + 'H',
			0,1,0,
			i*xres, 0, ptl,
			i*xres, 0, bz,
			(i+1)*xres, 0, bz, 0)
		t2 = struct.pack('<' + 'f'*12 + 'H',
			0,1,0,
			(i+1)*xres, 0, ptr,
			i*xres, 0, ptl,
			(i+1)*xres, 0, bz, 0)
		b1 = struct.pack('<' + 'f'*12 + 'H',
			0,-1,0,
			i*xres, (y-1)*yres, bz,
			i*xres, (y-1)*yres, pbl,
			(i+1)*xres, (y-1)*yres, bz, 0)
		b2 = struct.pack('<' + 'f'*12 + 'H',
			0,-1,0,
			i*xres, (y-1)*yres, pbl,
			(i+1)*xres, (y-1)*yres, pbr,
			(i+1)*xres, (y-1)*yres, bz, 0)
		outfl.write(t1)
		outfl.write(t2)
		outfl.write(b1)
		outfl.write(b2)

	# Add the other two sides
	for j in range(y-1):
		ptl = dem_data[j,0]
		ptr = dem_data[j+1,0]
		pbl = dem_data[j,x-1]
		pbr = dem_data[j,x-1]
		l1 = struct.pack('<' + 'f'*12 + 'H',
			-1,0,0,
			0, j*yres, bz,
			0, j*yres, ptl,
			0, (j+1)*yres, bz, 0)
		l2 = struct.pack('<' + 'f'*12 + 'H',
			-1,0,0,
			0, j*yres, ptl,
			0, (j+1)*yres, ptr,
			0, (j+1)*yres, bz, 0)
		r1 = struct.pack('<' + 'f'*12 + 'H',
			1,0,0,
			(x-1)*xres, j*yres, pbl,
			(x-1)*xres, j*yres, bz,
			(x-1)*xres, (j+1)*yres, bz, 0)
		r2 = struct.pack('<' + 'f'*12 + 'H',
			1,0,0,
			(x-1)*xres, (j+1)*yres, pbl,
			(x-1)*xres, j*yres, pbr,
			(x-1)*xres, (j+1)*yres, bz, 0)
		outfl.write(r1)
		outfl.write(r2)
		outfl.write(l1)
		outfl.write(l2)

	# Add a bottom
	f1 = struct.pack('<' + 'f'*12 + 'H',
			0,0,-1,
			0, 0, bz,
			0, (y-1)*yres, bz,
			(x-1)*xres, 0, bz, 0)
	f2 = struct.pack('<' + 'f'*12 + 'H',
			0,0,-1,
			0, (y-1)*yres, bz,
			(x-1)*xres, (y-1)*yres, bz,
			(x-1)*xres, 0, bz, 0)
	outfl.write(f1)
	outfl.write(f2)
		
outfl.close()
