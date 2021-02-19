import numpy as np
import argparse, os, struct
import rasterio as rio
import scipy.signal as spsig

def genFacets(grid, dx, dy,):
	h, w = grid.shape
	# Elev of bottom - min minus 10pct of max-min
	bz = grid.min() - (0.1 * (grid.max() - grid.min()))-1
	
	#nfacet = (w-1)*(h-1)*2
	# Number of facets for making a volume
	nfacet = (
		2 * (w - 1) * (h - 1) + 4 * (w - 1) + 4 * (h - 1) + 2
	)

	# Array to hold facet data
	# Cols 1-9 hold (x1,y1,z1,x2,...,z3) of the facet corners
	# Col 10-12 hold the normal vector to the facet
	f = np.zeros((nfacet, 12))

	# Build top surfce
	c = 0
	for i in range(h-1):
		for j in range(w-1):
			## Facet 1 - upper left, bottom left, bottom right
			# Point 1
			f[c,3] = j*dx
			f[c,4] = i*dy
			f[c,5] = grid[i,j]
			# Point 2
			f[c,6] = j*dx
			f[c,7] = (i+1)*dy
			f[c,8] = grid[i+1,j]
			# Point 3
			f[c,9] = (j+1)*dx
			f[c,10] = (i+1)*dy
			f[c,11] = grid[i+1,j+1]
			# Normal (p2-p1 X p3-p2)
			f[c,0:3] = -1*np.cross(f[c,6:9] - f[c,3:6], f[c,9:12]-f[c,6:9])
			c += 1

			## Facet 2 - upper left, upper right, bottom right
			# Point 1
			f[c,3] = j*dx
			f[c,4] = i*dy
			f[c,5] = grid[i,j]
			# Point 2
			f[c,6] = (j+1)*dx
			f[c,7] = (i+1)*dy
			f[c,8] = grid[i+1,j+1]
			# Point 3
			f[c,9] = (j+1)*dx
			f[c,10] = i*dy
			f[c,11] = grid[i,j+1]
			# Normal (p2-p3 X p1-p2)
			f[c,0:3] = np.cross(f[c,6:9] - f[c,9:12], f[c,3:6]-f[c,6:9])
			c += 1

	## Build sides
	# Left/Right sides
	j = 0
	norm = [-1,0,0]
	for i in range(h-1):
		## Facet 1 - upper left, bottom left, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = j*dx
		f[c,7] = i*dy
		f[c,8] = bz
		# Point 3
		f[c,9] = j*dx
		f[c,10] = (i+1)*dy
		f[c,11] = bz
		# Normal
		f[c,0:3] = norm
		c += 1

		## Facet 2 - upper left, upper right, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 3
		f[c,6] = j*dx
		f[c,7] = (i+1)*dy
		f[c,8] = bz
		# Point 2
		f[c,9] = j*dx
		f[c,10] = (i+1)*dy
		f[c,11] = grid[i+1,j]
		# Normal (p2-p3 X p1-p2)
		f[c,0:3] = norm
		c += 1

	j = w-1
	norm = [1,0,0]
	for i in range(h-1):
		## Facet 1 - upper left, bottom left, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = j*dx
		f[c,7] = (i+1)*dy
		f[c,8] = bz
		# Point 3
		f[c,9] = j*dx
		f[c,10] = i*dy
		f[c,11] = bz
		# Normal
		f[c,0:3] = norm
		c += 1

		## Facet 2 - upper left, upper right, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = j*dx
		f[c,7] = (i+1)*dy
		f[c,8] = grid[i+1,j]
		# Point 3
		f[c,9] = j*dx
		f[c,10] = (i+1)*dy
		f[c,11] = bz
		# Normal (p2-p3 X p1-p2)
		f[c,0:3] = norm
		c += 1

	# Top/Bottom sides
	i = 0
	norm = [0,1,0]
	for j in range(w-1):
		## Facet 1 - upper left, bottom left, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = (j+1)*dx
		f[c,7] = i*dy
		f[c,8] = bz
		# Point 3
		f[c,9] = j*dx
		f[c,10] = i*dy
		f[c,11] = bz
		# Normal
		f[c,0:3] = norm
		c += 1

		## Facet 2 - upper left, upper right, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = (j+1)*dx
		f[c,7] = i*dy
		f[c,8] = grid[i,j+1]
		# Point 3
		f[c,9] = (j+1)*dx
		f[c,10] = i*dy
		f[c,11] = bz
		# Normal (p2-p3 X p1-p2)
		f[c,0:3] = norm
		c += 1

	i = h-1
	norm = [0,-1,0]
	for j in range(w-1):
		## Facet 1 - upper left, bottom left, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = j*dx
		f[c,7] = i*dy
		f[c,8] = bz
		# Point 3
		f[c,9] = (j+1)*dx
		f[c,10] = i*dy
		f[c,11] = bz
		# Normal
		f[c,0:3] = norm
		c += 1

		## Facet 2 - upper left, upper right, bottom right
		# Point 1
		f[c,3] = j*dx
		f[c,4] = i*dy
		f[c,5] = grid[i,j]
		# Point 2
		f[c,6] = (j+1)*dx
		f[c,7] = i*dy
		f[c,8] = bz
		# Point 3
		f[c,9] = (j+1)*dx
		f[c,10] = i*dy
		f[c,11] = grid[i,j+1]
		# Normal (p2-p3 X p1-p2)
		f[c,0:3] = norm
		c += 1

	# Add bottom
	f[c,:] = [0,0,-1,0,0,bz,(w-1)*dx,0,bz,(w-1)*dx,(h-1)*dy,bz]
	f[c+1,:] = [0,0,-1,0,0,bz,(w-1)*dx,(h-1)*dy,bz,0,(h-1)*dy,bz]

	return f

def buildSTL(grid, xres, yres, stlFile):
	f = genFacets(grid, xres, yres)

	# Start making STL file
	fd = open(stlFile, "wb")

	head = struct.pack(
		"<" + "d" * 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	)  # 80 byte padding at beginning
	fd.write(head)

	nfac = struct.pack("<I", f.shape[0])
	fd.write(nfac)

	for i in range(f.shape[0]):
		row = struct.pack('<'+'f'*12, f[i, 0],f[i, 1],f[i, 2],f[i, 3],f[i, 4],f[i, 5],f[i, 6],f[i, 7],f[i, 8],f[i, 9],f[i, 10],f[i, 11])
		bc = struct.pack('<H', 0)
		fd.write(row+bc)

	fd.close()

	return 0


def main():
	# Set up CLI
	parser = argparse.ArgumentParser(
		description="""Program for generating a STL file from a dem
					   requires a square or rectangular DEM with no
					   no-data swaths. Horizontal units must match
					   vertical units. Output filename is auto-generated
					   if not specified.""",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	)
	parser.add_argument("DEM", help="Digital elevation model")
	parser.add_argument("-o", "--output", help="Output STL file", default=None)
	parser.add_argument(
		"-s",
		"--scale",
		help="Map scale of output STL will be 1:SCALE. Assumes DEM coordinates are meters.",
		type=float,
		default=1,
	)
	parser.add_argument(
		"-vx",
		"--vertical-exaggeration",
		help="Vertical exaggeration of output STL relative to input DEM",
		type=float,
		default=1,
	)
	args = parser.parse_args()

	# Open DEM and load data/metadata
	ds = rio.open(args.DEM, "r")
	grid = ds.read(1)
	w = ds.width
	h = ds.height
	dx = ds.transform[0]
	dy = ds.transform[4]
	ds.close()

	# Generate output file name if not specified
	if args.output == None:
		fext = os.path.splitext(args.DEM)[1]
		args.output = args.DEM.replace(
			fext,
			"_1to" + str(args.scale) + "_" + str(args.vertical_exaggeration) + "x.stl",
		)

	# Calculate scaled dimensions and resolution
	txdim = np.abs(dx * w)  # True dims in m
	tydim = np.abs(dy * h)
	sf = args.scale / 1000  # scale factor
	sxdim = txdim / sf  # Scaled dims in mm
	sydim = tydim / sf
	sdx = dx / sf  # Scale dx and dy
	sdy = dy / sf

	# Decimate points to about 1 per mm if necessary
	# Gets rid of extra data and noise from the stl
	if sdx < 1:
		xf = np.abs(int(np.ceil(1 / sdx)))
		grid = spsig.decimate(grid, xf, axis=1, zero_phase=True)
		sdx = sdx * xf
		sxdim = np.abs(sdx * grid.shape[1])
	if sdy < 1:
		yf = np.abs(int(np.ceil(1 / sdy)))
		grid = spsig.decimate(grid, yf, axis=0, zero_phase=True)
		sdy = sdy * yf
		sydim = np.abs(sdy * grid.shape[0])

	# Apply vx and scale factor to grid Z
	grid = grid * args.vertical_exaggeration/sf

	# Make stl file
	buildSTL(grid, sdx, sdy, args.output)

	return 0


main()
