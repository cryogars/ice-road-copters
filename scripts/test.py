import rasterio
from rasterio.merge import merge
from rasterio.plot import show
import glob
import os

# List all files in (to be variable)
search = glob.glob(os.path.join('./usgs_1m_tiles/MCS/*.tif'))

# Make a list of files
src_files_to_mosaic = []
for fp in search:
    src = rasterio.open(fp)
    src_files_to_mosaic.append(src)

# Merge list of tif files
mosaic, out_trans = merge(src_files_to_mosaic)

# Copy metadata to mosaic
out_meta = src.meta.copy()

# Output mosaic tif.. this will next be used as reference in pc_align
with rasterio.open('./out.tif', "w", **out_meta) as dest:
    dest.write(mosaic)

# Use ASP to convert from geoid to ellipsoid
dem_geoid('mosaic.tif --geoid NAVD88 --reverse-adjustment -o dem_ellipsoid')

