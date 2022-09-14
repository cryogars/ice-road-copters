# Import libraries
import os
import re
import sys
from os.path import exists, join, basename, dirname, abspath, isdir
import glob
from unittest import result
import geopandas as gpd
from laz2dem import iceroad_logging, cl_call
import json
import logging

log = logging.getLogger(__name__)

def clip_align(input_laz, buff_shp, result_dir, json_dir, log, dem_is_geoid, asp_dir, final_tif):
        # Clip clean_PC to the transform_area using PDAL
        # input_laz = join(result_dir, basename(in_dir)+'_unaligned.laz')
        clipped_pc = join(result_dir, 'clipped_PC.laz')
        json_path = join(json_dir, 'clip_to_shp.json')

        # Create .json file for PDAL clip
        json_pipeline = {
            "pipeline": [
                input_laz,
                {
                    "type":"filters.overlay",
                    "dimension":"Classification",
                    "datasource":buff_shp,
                    "layer":"buffered_area",
                    "column":"CLS"
                },
                {
                    "type":"filters.range",
                    "limits":"Classification[42:42]"
                },
                clipped_pc
            ]
        }
        with open(json_path,'w') as outfile:
            json.dump(json_pipeline, outfile, indent = 2)

        cl_call(f'pdal pipeline {json_path}', log)               

        # Check to see if output clipped point cloud was created
        if not exists(clipped_pc):
            log.info('Output point cloud not created')
            return -1

        log.info('Point cloud clipped to area')

        # Define paths for next if statement
        in_dem = join(result_dir, 'dem.tif')
        
        if dem_is_geoid is True:
            # ASP needs NAVD88 conversion to be in NAD83 (not WGS84)
            nad83_dem = join(result_dir, 'demNAD_tmp.tif')
            gdal_func = join(asp_dir, 'gdalwarp')
            cl_call(f'{gdal_func} -t_srs EPSG:26911 {in_dem} {nad83_dem}', log)
            # Use ASP to convert from geoid to ellipsoid
            ellisoid_dem = join(result_dir, 'dem_wgs')
            geoid_func = join(asp_dir, 'dem_geoid')
            cl_call(f'{geoid_func} --nodata_value -9999 {nad83_dem} \
                    --geoid NAVD88 --reverse-adjustment -o {ellisoid_dem}', log)
            # Set it back to WGS84
            ref_dem = join(result_dir, 'ellipsoid_DEM.tif')
            cl_call(f'{gdal_func} -t_srs EPSG:32611 {ellisoid_dem}-adj.tif {ref_dem}', log)

            # check for success
            if not exists(ref_dem):
                log.info('Conversion to ellipsoid failed')
                return -1

            log.info('Merged DEM converted to ellipsoid per user input')

        else:
            # cl_call('cp '+ in_dem +' '+ ref_dem, log)
            ref_dem = in_dem
            log.info('Merged DEM was kept in original ellipsoid form...')

        # Call ASP pc_align function on road and DEM and output translation/rotation matrix
        pc_align_func = join(asp_dir, 'pc_align')
        align_pc = join(result_dir,'pc-align',basename(final_tif))
        log.info('Beginning pc_align function...')
        cl_call(f'{pc_align_func} --max-displacement 5 --highest-accuracy \
                    {ref_dem} {clipped_pc} -o {align_pc}', log) # change run to something better
        
        # Apply transformation matrix to the entire laz and output points
        # https://groups.google.com/g/ames-stereo-pipeline-support/c/XVCJyXYXgIY/m/n8RRmGXJFQAJ
        transform_pc = join(result_dir,'pc-transform',basename(final_tif))
        cl_call(f'{pc_align_func} --max-displacement -1 --num-iterations 0 \
                    --initial-transform {align_pc}-transform.txt \
                    --save-transformed-source-points                            \
                    {ref_dem} {input_laz}   \
                    -o {transform_pc}', log)

        # Grid the output to a 0.5 meter tif (NOTE: this needs to be changed to 1m if using py3dep)
        point2dem_func = join(asp_dir, 'point2dem')
        # final_tif = join(ice_dir, 'pc-grid', 'run')
        cl_call(f'{point2dem_func} {transform_pc}-trans_source.laz \
                    --dem-spacing 0.5 --search-radius-factor 2 -o {final_tif}', log)
    
        return final_tif + '-DEM.tif'

# Find transformations/rotations via iceyroads and apply to whole point cloud
def laz_align(in_dir, 
            log,
            input_laz,
            canopy_laz,
            align_shp = 'transform_area/hwy_21/hwy_21_utm_edit_v2.shp',
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
    in_dir = abspath(in_dir)
    assert isdir(in_dir), 'in_dir must be directory'
    ice_dir = join(in_dir, 'ice-road')
    result_dir = join(ice_dir, 'results')
    json_dir = join(ice_dir,'jsons')

    # log = iceroad_logging(join(work_dir, 'logs'), debug = True, log_prefix='asp_align')
    log.info('Starting ASP align')

    # todo: since the buffer is in meters, need to ensure inputs are in UTM and same
    # Read in transform area (ice roads)
    log.info(f'Loading in shapefile {align_shp}')
    gdf = gpd.read_file(align_shp)

    # Buffer geom based on user input
    gdf['geometry'] = gdf.geometry.buffer(buffer_meters)

    # Create a new attribute to be used for PDAL clip/overlay
    gdf['CLS'] = 42

    # Save buffered shpfile to directory we just made
    buff_shp = join(result_dir, 'buffered_area.shp')
    gdf.to_file(buff_shp)

    # asp_dir = join(result_dir, 'asp')
    # os.makedirs(asp_dir, exist_ok= True)

    snow_final_tif = join(ice_dir, basename(in_dir)+'-snow')
    canopy_final_tif = join(ice_dir, basename(in_dir)+'-snow')
    if exists(snow_final_tif + '-DEM.tif'):
        while True:
            ans = input("Aligned tif already exists. Enter y to overwrite and n to use existing:")
            if ans.lower() == 'n':
                return snow_final_tif + '-DEM.tif', canopy_final_tif + '-DEM.tif'
            elif ans.lower() == 'y':
                break
        
    snow_tif = clip_align(input_laz=input_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = snow_final_tif)

    canopy_tif = clip_align(input_laz=canopy_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = canopy_final_tif)

    # For some reason this is returning 1 when a product IS created..
    if not exists(snow_tif):
       log.info('No final product created')
       return -1

    return snow_tif, canopy_tif

if __name__ == '__main__':
    laz_align('/Users/brent/Documents/MCS/mcs0407/results')