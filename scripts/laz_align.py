# Import libraries
import os
from os.path import exists, join, basename, abspath, isdir
from unittest import result
import geopandas as gpd
from laz2dem import cl_call, filter_dem
import json
import logging

log = logging.getLogger(__name__)

def clip_align(
        input_laz: str, buff_shp: str, result_dir: str, json_dir: str, log: logging.Logger,
        dem_is_geoid: bool, asp_dir: str, final_tif: str, dem_fp: str, is_canopy: bool = False,
        las_extra_byte_format: bool = False, use_dem_filter: bool = True
) -> str:
    
    # Clip clean_PC to the transform_area using PDAL
    # input_laz = join(result_dir, basename(in_dir)+'_unaligned.laz')
    clipped_pc = join(result_dir, 'clipped_PC.laz')
    json_path = join(json_dir, 'clip_to_shp.json')

    if not exists(input_laz):
        raise FileNotFoundError(f"Input point cloud not found: {input_laz}")
    if not exists(buff_shp):
        raise FileNotFoundError(f"Buffered shapefile missing: {buff_shp}")

    if not asp_dir or not exists(asp_dir):
        raise RuntimeError(f"ASP directory invalid or missing: {asp_dir}")


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

    
    if dem_is_geoid is True:
        # ASP needs NAVD88 conversion to be in NAD83 (not WGS84)
        nad83_dem = join(result_dir, 'demNAD_tmp.tif')
        gdal_func = join(asp_dir, 'gdalwarp')
        cl_call(f'{gdal_func} -t_srs EPSG:26911 {dem_fp} {nad83_dem}', log)
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
            raise Exception('Conversion to ellipsoid failed')

        log.info('Merged DEM converted to ellipsoid per user input')

        # Cleanup temporary files
        log.info("Cleaning up temporary files")
        for tmp in (nad83_dem, ellisoid_dem + "-adj.tif"):
            if exists(tmp):
                os.remove(tmp)
        log.info("Cleanup finished!")

    else:
        ref_dem = dem_fp
        log.info('Merged DEM was kept in original ellipsoid form...')

    # Call ASP pc_align function on road and DEM and output translation/rotation matrix
    align_pc = join(result_dir,'pc-align',basename(final_tif))
    pc_align_func = join(asp_dir, 'pc_align') 
    log.info('Beginning pc_align function...')
    cl_call(f'{pc_align_func} --max-displacement 5 --highest-accuracy \
                {ref_dem} {clipped_pc} -o {align_pc}', log)
    
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
        

    # Apply transformation matrix to the entire laz and output points
    # https://groups.google.com/g/ames-stereo-pipeline-support/c/XVCJyXYXgIY/m/n8RRmGXJFQAJ
    transform_pc = join(result_dir,'pc-transform',basename(final_tif))
    cl_call(f'{pc_align_func} --max-displacement -1 --num-iterations 0 \
                --initial-transform {align_pc}-transform.txt \
                --save-transformed-source-points                            \
                {ref_dem} {input_laz}   \
                -o {transform_pc}', log)

    # Apply DEM filter to aligned point cloud
    aligned_laz = f'{transform_pc}-trans_source.laz'
    if not exists(aligned_laz):
        raise Exception(f"Aligned point cloud not created: {aligned_laz}")
    
    if use_dem_filter and not is_canopy:
        # Run post-alignment DEM-based filtering via PDAL
        log.info("Post-alignment DEM filtering ENABLED for for ground/DTM")
        stem = basename(final_tif).replace(".tif", "")
        filtered_laz = join(result_dir, 'pc-transform', f"{stem}_filtered.laz")
        filtered_laz = filter_dem(
            input_laz=aligned_laz,
            dem_fp=ref_dem,
            outlas=filtered_laz,
            json_dir=json_dir,
            log=log,
        )
    else:
        log.info("Post-alignment DEM filtering not applied")
        filtered_laz = aligned_laz

    # Grid the output to a 0.5 meter tif (NOTE: this needs to be changed to 1m if using py3dep)
    point2dem_func = join(asp_dir, 'point2dem')
    cl_call(f'{point2dem_func} {filtered_laz} \
                --dem-spacing 0.5 --search-radius-factor 2 -o {final_tif}', log)

    return final_tif + '-DEM.tif'




# Find transformations/rotations via iceyroads and apply to whole point cloud
def laz_align(in_dir: str, 
            log: logging.Logger,
            input_laz: str,
            canopy_laz: str,
            user_dem: str,
            align_shp: str = 'transform_area/hwy_21/hwy_21_utm_edit_v2.shp',
            buffer_meters: float = 3.0, 
            dem_is_geoid: bool = False, 
            asp_dir: str | None = None,
            las_extra_byte_format: bool=False,
            use_dem_filter : bool =True,
            force_overwrite: bool = False
            ) -> tuple[str, str]:
    """
    Align point cloud using snow-off road polygon.

    Parameters:
    in_dir: filepath to directory to run in
    align_shp: filepath to shapefile to clip point cloud to
    buffer_meters: number of meters to buffer geometry
    dem_is_geoid: leave as geoid or convert to ellispoid
    asp_dir: filepath to ASP bin directory
    use_dem_filter: apply post-alignment DEM filtering
    user_dem: user-specified DEM (required)
    force_overwrite: force overwrite aligned TIF (default = False)

    Returns:
    final_tif: filepath to output corrected point cloud
    """
    in_dir = abspath(in_dir)
    assert isdir(in_dir), 'in_dir must be directory'
    ice_dir = join(in_dir, 'ice-road')
    result_dir = join(ice_dir, 'results')
    json_dir = join(ice_dir,'jsons')

    os.makedirs(ice_dir, exist_ok= True)
    os.makedirs(result_dir, exist_ok= True)
    os.makedirs(json_dir, exist_ok= True)
    
    #------------
    # Process DEM
    #------------

    if not user_dem.lower().endswith((".tif", ".tiff")):
        raise RuntimeError(f"User DEM must be a .tif file: {user_dem}")

    dem_fp = join(result_dir, "dem.tif")
    log.info(f"Using user provided DEM: {user_dem}")
    cl_call(f"cp {user_dem} {dem_fp}", log)

    if not exists(dem_fp):
        raise RuntimeError(f"Failed to copy DEM to destination: expected file at {dem_fp}")
    
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
    if exists(snow_final_tif + '.tif') and exists(canopy_final_tif + '.tif') and not force_overwrite:
         log.info("Reusing existing aligned TIF")
         return snow_final_tif + '.tif', canopy_final_tif+ '.tif'
    elif exists(snow_final_tif + '.tif') and exists(canopy_final_tif + '.tif') and force_overwrite:
        log.info("Overwriting aligned TIF")
        
    snow_tif = clip_align(input_laz=input_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = snow_final_tif, is_canopy=False, las_extra_byte_format=las_extra_byte_format,
        dem_fp=dem_fp, use_dem_filter=use_dem_filter)

    canopy_tif = clip_align(input_laz=canopy_laz, buff_shp=buff_shp, result_dir=result_dir,\
        json_dir=json_dir, log = log, dem_is_geoid=dem_is_geoid, asp_dir=asp_dir,\
        final_tif = canopy_final_tif, is_canopy=True, las_extra_byte_format=las_extra_byte_format,
        dem_fp=dem_fp, use_dem_filter=False)

    # For some reason this is returning 1 when a product IS created..
    if not exists(snow_tif):
       log.info(f'Can not find {snow_tif}')
       raise Exception('No final product created')

    return snow_tif, canopy_tif
