# To run For now just run this script with no params

# Import libraries
import os

import glob
import geopandas as gpd
import whitebox
from laz2dem import iceroad_logging

# Find transformations/rotations via iceyroads and apply to whole point cloud
def laz_align(transform_area='hwy_21', buffer_meters=2.5, geoid=False):
    '''
    to be ran after Zach's code... transform_area = 'hwy_21' for now
    '''
    log = iceroad_logging('./data/logs', debug = True, log_prefix='laz_align')
    # Hard code in /data/results as the directory
    # This works as long as user supplies las/laz in `data`... all following zach's code
    dirname =  os.path.abspath('./data/results/')

    # TODO: since the buffer is in meters, need to ensure inputs are in UTM and same
    path = './transform_area/'+transform_area+'/*.shp'
    for filename in glob.glob(path):
        # Read in transform area (ice roads)
        gdf = gpd.read_file(filename)
        # Buffer geom based on user input
        gdf['geometry'] = gdf.geometry.buffer(buffer_meters)
        # Save buffered shpfile to directory we just made
        gdf.to_file(dirname+'/buffered_area.shp')

    # Clip clean_PC to the transform_area using whitebox-python
    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = dirname
    wbt.clip_lidar_to_polygon(i='data.laz', 
                              polygons='buffered_area.shp',
                              output='clipped_PC.laz')
    log.info('PC successfully clipped to area')

    if not geoid:
        # ASP needs NAVD88 conversion to be in NAD83 (not WGS84)
        os.system('./ASP/bin/gdalwarp -t_srs EPSG:26911 ./data/results/dem.tif ./data/results/demNAD_tmp.tif')

        # Use ASP to convert from geoid to ellipsoid
        os.system('./ASP/bin/dem_geoid --nodata_value -9999 ./data/results/demNAD_tmp.tif \
                   --geoid NAVD88 --reverse-adjustment -o ./data/results/dem_wgs')

        # Set it back to WGS84
        os.system('./ASP/bin/gdalwarp -t_srs EPSG:32611 ./data/results/dem_wgs-adj.tif ./data/results/ref_PC.tif')
        log.info('Merged DEM converted to ellipsoid per user input')

    else:
        log.info('Merged DEM was kept as geoid per user input')


    # Call ASP pc_align function on road and DEM and output translation/rotation matrix
    log.info('Beginning pc_align function...')
    os.system('./ASP/bin/pc_align --max-displacement 5 --highest-accuracy \
                ./data/results/ref_PC.tif ./data/results/clipped_PC.laz -o ./data/results/pc-align/run')
    
    # Apply transformation matrix to the entire laz and output points
    # https://groups.google.com/g/ames-stereo-pipeline-support/c/XVCJyXYXgIY/m/n8RRmGXJFQAJ
    os.system('./ASP/bin/pc_align --max-displacement -1 --num-iterations 0 \
                --initial-transform data/results/pc-align/run-transform.txt \
                --save-transformed-source-points                            \
                ./data/results/ref_PC.tif ./data/results/data.laz   \
                -o ./data/results/pc-transform/run')

    # Grid the output to a 1 meter tif
    os.system('./ASP/bin/point2dem ./data/results/pc-transform/run-trans_source.laz \
                --dem-spacing 1 --search-radius-factor 2 -o ./data/results/pc-grid/run')


# To run For now just run this script
laz_align()