# Import libraries
import os

import shutil
import glob
from datetime import datetime
import geopandas as gpd
import whitebox
import rasterio
from rasterio.merge import merge
from rasterio.plot import show

# This is the main
def iceroads_asp(clean_PC, transform_area, geoid=False, buffer_meters=5):
    '''
    ***********INPUTS***********

    `clean_PC` = name of lidar file with cleaned, ground classified point cloud (.laz). Should be saved in data_processed

    `transform_area` = folder name containing small subset area used to align `clean_PC` with the USGS 3DEP 1m DEM (.shp)

    `geoid` = Does your input lidar point cloud have a vertical datum that is geoid? (T/F). Default is False.

    `buffer_meters` = the linear buffer [meters] you would like to apply to `transform_area`. Default is 5 meters.


    ***********OUTPUTS***********

    All output files will be stored under data_processed in the time function is called.

    `clipped_PC` = This is the first output generated from the function and contains only points within the user specified buffer.
    
    `aligned_PC` = xxx

    `alignment_stats` = xxx

    `alignment matrix` = xxx

    other ASP outputs like kmz file virtual gcp...


    '''

    # Create processing directory with datetime.now
    abs = os.path.abspath('./data_processed/')
    today = datetime.now()
    dirname = abs+'/iceroads_asp_' + today.strftime('%Y%m%d%H%M')+ '_' + str(buffer_meters)+'m_Buffer'

    # Delete if already exists (happens if run more than once in a single minute)
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)

    # Create directory on disk
    os.mkdir(dirname)

    # Right now this only supports one shapefile
    # But this one shapefile can have many elements to it..
    # But potentially we could build in functionality here to combine shps here..
    path = './transform_area/'+transform_area+'/*.shp'
    for filename in glob.glob(path):
        
        # Read in transform area (ice roads)
        gdf = gpd.read_file(filename)
        
        # Buffer geom based on user input
        gdf['geometry'] = gdf.geometry.buffer(buffer_meters)

        # Save buffered shpfile to directory we just made
        gdf.to_file(dirname+'/buffered_area.shp')

    # Create a symlink for whitebox tools
    src = abs + '/' + clean_PC
    dst = dirname+'/clean_PC.laz'
    os.symlink(src, dst)
  
    # Clip clean_PC to the transform_area using whitebox-python
    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = dirname
    wbt.clip_lidar_to_polygon(i='clean_PC.laz', 
                              polygons='buffered_area.shp',
                              output='clipped_PC.laz')
    print('[INFO] PC successfully clipped to area')

    # Load in all of the tiles and prep for ASP....
    print('[INFO] USGS TIF files being loaded and merged...')

    # TODO: COPY AND PASTE TEST HERE

    # Use ASP to convert from geoid to ellipsoid
    if not geoid:
        #dem_geoid('mosaic.tif --geoid NAVD88 --reverse-adjustment -o dem_ellipsoid')
        print('[INFO] Merged DEM converted to ellipsoid per user input')
    else:
        print('[INFO] Merged DEM was kept as geoid per user input')


    # Call ASP pc_align function on road and DEM
    print('[INFO] Beginning pc_align function...')
    #pc_align()
    
    
    # Apply transformation matrix to laz file..Calculate RMSE based on changes...



    # Export laz as aligned_PC

# Used for testing
iceroads_asp(clean_pc='clean_PC.laz', transform_area='hwy_21', geoid=False, buffer_meters=5)