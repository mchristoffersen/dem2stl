# dem2stl
dem2stl is a python script to transform a square or rectangular DEM into a 3D print-able STL file. Argparse printout explaining usage - 

usage: dem2stl.py [-h] [-o OUTPUT] [-s SCALE] [-vx VERTICAL_EXAGGERATION]
                  [-p POINT_DENSITY]
                  DEM

Program for generating a STL file from a dem requires a square or rectangular
DEM with no no-data swaths. Horizontal units must match vertical units. Output
filename is auto-generated if not specified.

positional arguments:
  DEM                   Digital elevation model

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output STL file (default: None)
  -s SCALE, --scale SCALE
                        Map scale of output STL will be 1:SCALE. Assumes DEM
                        coordinates are meters. (default: 1)
  -vx VERTICAL_EXAGGERATION, --vertical-exaggeration VERTICAL_EXAGGERATION
                        Vertical exaggeration of output STL relative to input
                        DEM (default: 1)
  -p POINT_DENSITY, --point-density POINT_DENSITY
                        Output STL point density in points per mm (default: 1)
