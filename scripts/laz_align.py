# Import libraries
import os

import glob
import geopandas as gpd
import whitebox

# Find transformations/rotations via iceyroads and apply to whole point cloud
def laz_align(transform_area='hwy_21', buffer_meters=5, geoid=False):
    '''
    to be ran after Zach's code... transform_area = 'hwy_21' for now
    '''
    
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
    wbt.clip_lidar_to_polygon(i='merge.laz', 
                              polygons='buffered_area.shp',
                              output='clipped_PC.laz')
    print('[INFO] PC successfully clipped to area')

    if not geoid:
        # ASP needs NAVD88 conversion to be in NAD83 (not WGS84)
        os.system('./ASP/bin/gdalwarp -t_srs EPSG:26911 ./data/results/dem.tif ./data/results/demNAD_tmp.tif')

        # Use ASP to convert from geoid to ellipsoid
        os.system('./ASP/bin/dem_geoid --nodata_value -9999 ./data/results/demNAD_tmp.tif \
                   --geoid NAVD88 --reverse-adjustment -o ./data/results/dem_wgs')

        # Set it back to WGS84
        os.system('./ASP/bin/gdalwarp -t_srs EPSG:32611 ./data/results/dem_wgs-adj.tif ./data/results/ref_PC.tif')
        print('[INFO] Merged DEM converted to ellipsoid per user input')

    else:
        print('[INFO] Merged DEM was kept as geoid per user input')


    # Call ASP pc_align function on road and DEM and output translation/rotation matrix
    print('[INFO] Beginning pc_align function...')
    os.system('./ASP/bin/pc_align --max-displacement 5 --highest-accuracy \
                --save-transformed-source-points --save-inv-transformed-reference-points \
                ./data/results/ref_PC.tif ./data/results/clipped_PC.laz -o ./data/results/pc-align/run')
    
    
    # Apply transformation matrix to laz file (FINAL)..Calculate RMSE based on changes...


# For Testing
laz_align()