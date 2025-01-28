# Import libraries
from os.path import exists, join, basename, dirname, abspath, isdir
from unittest import result
import geopandas as gpd
from laz2dem import cl_call
import json
import logging

log = logging.getLogger(__name__)

def clip_align(input_laz, buff_shp, result_dir, json_dir, log, dem_is_geoid, asp_dir, final_tif, is_canopy=False, las_extra_byte_format=False):
    
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
        raise Exception('Output point cloud not created')

    log.info('Point cloud clipped to area')

    # Define paths
    ref_dem = join(result_dir, 'dem.tif')

    # Call ASP pc_align function on road and DEM and output translation/rotation matrix
    align_pc = join(result_dir,'pc-align',basename(final_tif))
    pc_align_func = join(asp_dir, 'pc_align') 
    log.info('Beginning pc_align function...')
    cl_call(f'{pc_align_func} --max-displacement 5 --highest-accuracy \
                {ref_dem} {clipped_pc} -o {align_pc}', log)
    
    # FOR GRAIN SIZE CALCS
    # Since there are issues in transforming the point cloud and retaining reflectance,
    # the best I can do is translation only and no rotation..
    # Therefore, in this section, if the mode is set to calc Grain Size, an additional pc_align 
    # will be called in order to save the X,Y,Z translation only. This will not be applied
    # to the snow depth products, so there may be some subtle differences when comparing between the two. 
    # However, this is in hopes to retain the higher information where we can..
    # --compute-translation-only
    if las_extra_byte_format is True and is_canopy is False:
        transform_pc_temp = join(result_dir,'pc-align-translation-only','temp')
        cl_call(f'{pc_align_func} --max-displacement 5 --highest-accuracy \
                --compute-translation-only   \
                    {ref_dem} {clipped_pc}   \
                    -o {transform_pc_temp}', log)     
        
    # FOR NORMAL ICE-ROAD PIPELINE (includes rotation)
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
            buffer_meters=3.0, 
            dem_is_geoid=False, 
            asp_dir = None,
            las_extra_byte_format=False):
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

    # Buffer geom based on user input. NOTE: we assume buffer_meters is the entire width. 
    # So, must divide by 2 here to get the right distance from centerline.
    log.info(f'Buffer width of {buffer_meters} m is used. This is {buffer_meters / 2} m from centerline.')
    gdf['geometry'] = gdf.geometry.buffer(buffer_meters / 2)

    # Create a new attribute to be used for PDAL clip/overlay
    gdf['CLS'] = 42

    # Save buffered shpfile to directory we just made
    buff_shp = join(result_dir, 'buffered_area.shp')
    gdf.to_file(buff_shp)

    # asp_dir = join(result_dir, 'asp')
    # os.makedirs(asp_dir, exist_ok= True)

    snow_final_tif = join(ice_dir, basename(in_dir)+'-snow')
    canopy_final_tif = join(ice_dir, basename(in_dir)+'-canopy')
    if exists(snow_final_tif + '.tif') and exists(canopy_final_tif + '.tif'):
        while True:
            ans = input("Aligned tif already exists. Enter y to overwrite and n to use existing:")
            if ans.lower() == 'n':
                return snow_final_tif + '.tif', canopy_final_tif+ '.tif'
            elif ans.lower() == 'y':
                break
        
    snow_tif = clip_align(input_laz=input_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = snow_final_tif, is_canopy=False, las_extra_byte_format=las_extra_byte_format)

    canopy_tif = clip_align(input_laz=canopy_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = canopy_final_tif, is_canopy=True, las_extra_byte_format=las_extra_byte_format)

    # For some reason this is returning 1 when a product IS created..
    if not exists(snow_tif):
       log.info(f'Can not find {snow_tif}')
       raise Exception('No final product created')

    return snow_tif, canopy_tif



if __name__ == '__main__':
    laz_align('/Users/brent/Documents/MCS/mcs0407/results')