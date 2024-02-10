# dem2stl
dem2stl is a python script to transform a square or rectangular DEM into a 3D print-able STL file. It uses (rasterio)[https://rasterio.readthedocs.io/en/stable/] to open DEMs, so should be able to handle any file format that rasterio can deal with.
  
dem2stl can be installed with pip:  
`pip install git+https://github.com/mchristoffersen/dem2stl.git`  

Once installed you should have `dem2stl` available at the command line.  
  
dem2stl can't handle DEMs with no-data cells. Also make sure that the elevation units are the same as the map projection units. Mixed unit DEMs (e.g. x/y in decimal degrees and z in meters) won't work well. The program reads in the first band of the raster it is fed, so if you have a multi-band file you'll either need to make sure that the elevation data you want to make a model from is the first band, modify dem2stl to read in the bad you need (either change the ds.read(1) line in main() or maybe add a command line option?), or generate a single band file from your multi band file and feed that to dem2stl.
  
Argparse printout explaining usage -  
```
usage: dem2stl [-h] [-o OUTPUT] [-s SCALE] [-vx VERTICAL_EXAGGERATION] [-p POINT_DENSITY] DEM

Program for generating a STL file from a dem requires a square or rectangular DEM with no no-data swaths.
Horizontal units must match vertical units. Output filename is auto-generated if not specified.

positional arguments:
  DEM                   Digital elevation model

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output STL file (default: None)
  -s SCALE, --scale SCALE
                        Map scale of output STL will be 1:SCALE. Assumes DEM coordinates are meters. (default:
                        1)
  -vx VERTICAL_EXAGGERATION, --vertical-exaggeration VERTICAL_EXAGGERATION
                        Vertical exaggeration of output STL relative to input DEM (default: 1)
  -p POINT_DENSITY, --point-density POINT_DENSITY
                        Output STL point density in points per mm (default: 1)
```
