# Test that the conda environment is working properly
import open3d
import pdal
from osgeo import gdal
import numpy
import pandas
import matplotlib
import geopandas
import descartes
import shapely
import contextily
import scipy
import haversine
import utm
import bokeh
import panel
import holoviews
import geoviews
import seaborn
import rasterstats
import hvplot
import rasterio
import xarray
import s3fs
import laspy

print('Conda env works!')

# Test that ASP is working properly
import os
os.system('./ASP/bin/pc_align --help')