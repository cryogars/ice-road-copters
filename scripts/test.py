import rasterio
from rasterio.merge import merge
from rasterio.plot import show
import glob
import os

# Note to self... still testing... having trouble downloading 1m dems from USGS
# still need to work out projection issue and finalize dem_geoid

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
os.system('./ASP/bin/dem_geoid out.tif --geoid NAVD88 --reverse-adjustment -o dem_ellipsoid.tif')

