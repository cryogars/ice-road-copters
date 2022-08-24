# To run For now just run this script with no params

# Import libraries
import os
import sys
from os.path import exists, join, basename, dirname, abspath, isdir
import glob
import geopandas as gpd
import whitebox
from laz2dem import iceroad_logging, cl_call
import logging

log = logging.getLogger(__name__)

# Find transformations/rotations via iceyroads and apply to whole point cloud
def laz_align(work_dir, 
            hwy_21_shp = '/Users/brent/Code/ice-road-copters/transform_area/hwy_21/hwy_21_utm_edit_v2.shp',
            buffer_meters=2.5, 
            dem_is_geoid=False, 
            asp_dir = None):
    '''
    Align point cloud using snow-off road polygon.

    Parameters:
    work_dir (str): filepath to directory to run in
    hwy_21_shp (str): filepath to shapefile to clip point cloud to
    buffer_meters (float): number of meters to buffer geometry
    geoid (bool): leave as geoid or convert to ellispoid
    asp_dir (str): filepath to ASP bin directory

    Returns:
    final_tif (str): filepath to output corrected point cloud
    '''
    if not asp_dir:
        asp_dir = abspath(join(dirname(dirname(work_dir)), 'ASP', 'bin'))
    work_dir = abspath(work_dir)
    assert isdir(work_dir), 'work_dir must be directory'

    # log = iceroad_logging(join(work_dir, 'logs'), debug = True, log_prefix='asp_align')
    log.info('Starting ASP align')

    # todo: since the buffer is in meters, need to ensure inputs are in UTM and same
    # Read in transform area (ice roads)
    gdf = gpd.read_file(hwy_21_shp)

    # Buffer geom based on user input
    gdf['geometry'] = gdf.geometry.buffer(buffer_meters)

    # Save buffered shpfile to directory we just made
    buff_shp = join(work_dir, 'buffered_area.shp')
    gdf.to_file(buff_shp)

    # Clip clean_PC to the transform_area using whitebox-python
    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = work_dir
    input_laz = join(work_dir, basename(dirname(work_dir))+'.laz')
    clipped_pc = join(work_dir, 'clipped_PC.laz')

    # If does exist give user option to clip lidar
    if exists(clipped_pc): 
        done = False
        while not done:
            ans = input("Clipped point cloud already exists. Enter y to overwrite and n to use existing:")
            if ans.lower() == 'n':
                done = True
            elif ans.lower() == 'y':
                wbt.clip_lidar_to_polygon(i=input_laz, 
                            polygons=buff_shp,
                            output=clipped_pc)
                done = True
    
    # If does not exist - make clipped lidar
    if not exists(clipped_pc):
        wbt.clip_lidar_to_polygon(i=input_laz, 
                                  polygons=buff_shp,
                                  output=clipped_pc)
    log.info('Point cloud clipped to area')

    # Define paths for next if statement
    in_dem = join(work_dir, 'dem.tif')
    ref_dem = join(work_dir, 'ref_PC.tif')
    
    if dem_is_geoid is True:
        # ASP needs NAVD88 conversion to be in NAD83 (not WGS84)
        nad83_dem = join(work_dir, 'demNAD_tmp.tif')
        gdal_func = join(asp_dir, 'gdalwarp')
        cl_call(f'{gdal_func} -t_srs EPSG:26911 {in_dem} {nad83_dem}', log)
        # Use ASP to convert from geoid to ellipsoid
        ellisoid_dem = join(work_dir, 'dem_wgs')
        geoid_func = join(asp_dir, 'dem_geoid')
        cl_call(f'{geoid_func} --nodata_value -9999 {nad83_dem} \
                   --geoid NAVD88 --reverse-adjustment -o {ellisoid_dem}', log)
        # Set it back to WGS84
        cl_call(f'{gdal_func} -t_srs EPSG:32611 {ellisoid_dem}-adj.tif {ref_dem}', log)

        # check for success
        if not exists(ref_dem):
            log.info('Conversion to ellipsoid failed')
            return 1

        log.info('Merged DEM converted to ellipsoid per user input')

    else:
        cl_call('cp '+ in_dem +' '+ ref_dem, log)
        log.info('Merged DEM was kept in original ellipsoid form...')

    # Call ASP pc_align function on road and DEM and output translation/rotation matrix
    pc_align_func = join(asp_dir, 'pc_align')
    align_pc = join(work_dir,'pc-align','run')
    log.info('Beginning pc_align function...')
    cl_call(f'{pc_align_func} --max-displacement 5 --highest-accuracy \
                {ref_dem} {clipped_pc} -o {align_pc}', log) # change run to somwthing better
    
    # Apply transformation matrix to the entire laz and output points
    # https://groups.google.com/g/ames-stereo-pipeline-support/c/XVCJyXYXgIY/m/n8RRmGXJFQAJ
    transform_pc = join(work_dir,'pc-transform','run')
    cl_call(f'{pc_align_func} --max-displacement -1 --num-iterations 0 \
                --initial-transform {align_pc}-transform.txt \
                --save-transformed-source-points                            \
                {ref_dem} {input_laz}   \
                -o {transform_pc}', log)

    # Grid the output to a 1 meter tif
    point2dem_func = join(asp_dir, 'point2dem')
    final_tif = join(work_dir, 'pc-grid', 'run')
    cl_call(f'{point2dem_func} {transform_pc}-trans_source.laz \
                --dem-spacing 0.5 --search-radius-factor 2 -o {final_tif}', log)

    if not exists(final_tif):
        log.info('No final product created')
        return 1

    return final_tif


# To run For now just run this script
if __name__ == '__main__':
    laz_align('/Users/brent/Documents/MCS/mcs0407/results')
