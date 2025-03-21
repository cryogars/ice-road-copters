#!/usr/bin/env python
"""
Takes input reference dem data (point cloud or raster) and will transform horizontal and vertical references.

Essentially a wrapper for ASP's dem_geoid but also allow for point clouds via PDAL.

Usage:
    transform_data.py <reference_data> [-e epsg] [-t transform_command] [-g geoid_model] [-a asp_dir] [-d debug]

Options:
    -e epsg                   Desired EPSG code you want (for example, 32611)
    -t transform_command      Either "to_geoid" or "to_ellipsoid".
    -g geoid_model            Either "EGM2008" or "NAVD88". Use EGM2008 if reference data is in WGS. Use NAVD88 if reference data in NAD.
    -a asp_dir                Directory with ASP binary files
    -d debug                  Turns on debugging logging

"""
import json
from docopt import docopt
from glob import glob
from os.path import abspath,exists, join, basename, dirname
import rasterio
import laspy
from pyproj import CRS
from datetime import datetime
import logging
import sys
import os

# local imports
from dir_space_strip import replace_white_spaces
from laz2dem import cl_call


def transform_pc(in_dir, epsg, asp_dir, geoid_dir, geoid_model, log):
    '''
    If user reference data is a point cloud, performs geoid transformation using PDAL.

    Inputs:
        in_dir - absolute path to point cloud file input by user.
        epsg - integer for EPSG code input by user.
        asp_dir - absolute path to ASP bin directory.

    Outputs:
        corrected_data - absolute path to transformed data.

    '''

    with laspy.open(in_dir) as las:
        hdr = las.header
        pc_crs = hdr.parse_crs()
        datum_ref = pc_crs.datum.name
    
    # Log info for user
    log.info(f'CRS for input data: {pc_crs}')
    log.info(f'~~~~~~~~~~~~~~~~~~~~')
    log.info(f'Input data vertical datum: {datum_ref}')

    # Ensure vertical datum of user data is not None
    if datum_ref is None:
        raise Exception(f'Vertical datum was found as {datum_ref}. Check log file to inspect CRS.')

    # geoid file to use from ASP
    if geoid_model.lower() == 'navd88':
        geoid_data_str = 'navd88.tif'
    elif geoid_model.lower() == 'egm2008':
        geoid_data_str = 'egm2008.jp2'
    else:
        raise Exception("Needs to be egm2008 or navd88.")

    # Path to the geoid from Ames Stereo
    parent_dir_asp = dirname(asp_dir)
    asp_geoid_data = join(parent_dir_asp, 'share', 'geoids', geoid_data_str)

    # corrected path to point cloud
    corrected_data = join(geoid_dir, f'reference-PC.{extension}')

    # Next, Make a PDAL json
    json_template = [
        f"{in_dir}",
        {
            "type":"filters.reprojection",
            "out_srs":f"+init=EPSG:{epsg} +geoidgrids={asp_geoid_data}"  
        },
        {
            "type":"writers.las",
            "a_srs":f"EPSG:{epsg}",
            "filename":f"{corrected_data}"
        }
    ]

    # Save json file
    json_to_use = join(geoid_dir, 'geoid_transformation.json')
    with open(json_to_use,'w') as outfile:
        json.dump(json_template, outfile, indent = 2)

    # Setting up command
    pipeline_cmd = f'pdal pipeline -i {json_to_use} -v 8'
    log.debug(f"Using pipeline command: {pipeline_cmd}")

    # run command
    cl_call(pipeline_cmd, log)


    return corrected_data


def transform_raster(in_dir, epsg, asp_dir, geoid_dir, geoid_model, transform_command, log ):
    '''
    If user reference data is a raster, performs geoid transformation using Ames Stereo Pipeline.

    Inputs:
        in_dir - absolute path to raster file input by user.
        epsg - integer for EPSG code input by user.
        asp_dir - absolute path to ASP bin directory.
        transform_command - either "to_geoid" or "to_ellipsoid".

    Outputs:
        corrected_data - absolute path to transformed data.

    '''
     # Check transform command and code to ASP commands
    if transform_command == 'to_ellipsoid':
        transform_cmd_for_asp = '--reverse-adjustment'
    elif transform_command == 'to_geoid':
        transform_cmd_for_asp = ''
    else:
        raise Exception(f'transform command must be either "to_geoid" or "to_ellipsoid"')

    # Load CRS and nodata from DEM
    ref_dem = rasterio.open(in_dir)
    nodata_value = ref_dem.nodata
    ref_dem_crs = ref_dem.crs.to_epsg()

    epsg_object = CRS.from_epsg(ref_dem_crs)
    datum_ref = epsg_object.datum.name

    # Log info for user
    log.info(f'CRS for input data: {ref_dem_crs}')
    log.info(f'~~~~~~~~~~~~~~~~~~~~')
    log.info(f'Input data vertical datum: {datum_ref}')

    # Ensure vertical datum of user data is not None
    if datum_ref is None:
        raise Exception(f'Vertical datum was found as {datum_ref}. Check log file to inspect CRS.')

    # Use ASP to convert
    # set up post tranform path
    transform_dem = join(geoid_dir, 'temp-transform')

    # set up calls to ASP
    geoid_func = join(asp_dir, 'dem_geoid')
    gdal_func = join(asp_dir, 'gdalwarp')

    # geoid file to use from ASP
    if geoid_model.lower() == 'navd88':
        geoid_cmd = 'NAVD88'
    elif geoid_model.lower() == 'egm2008':
        geoid_cmd = 'EGM2008'
    else:
        raise Exception("Needs to be egm2008 or navd88.")
    
    # Call geoid_func (dem_geoid)
    # NOTE: GEOID -> ELLIP === --reverse-adjustment ... other way is blank
    cl_call(f'{geoid_func} --nodata_value {nodata_value} {in_dir} \
            --geoid {geoid_cmd} {transform_cmd_for_asp} -o {transform_dem}', log)
    
    # Match CRS now to lidar (assign, it gets lost in ASP)
    corrected_data = join(geoid_dir, 'reference-DEM.tif')
    cl_call(f'{gdal_func} -t_srs EPSG:{epsg} {transform_dem}-adj.tif {corrected_data}', log)


    return corrected_data


if __name__ == '__main__':

    start_time = datetime.now()
    
    # get command line args
    args = docopt(__doc__)
    in_dir = args.get('<reference_data>')
    in_dir = abspath(in_dir)
    epsg = int(args.get('-e'))
    transform_command = args.get('-t')
    debug = args.get('-d')
    geoid_model = args.get('-g')
    asp_dir = args.get('-a')
    if asp_dir:
        asp_dir = abspath(asp_dir)
        if basename(asp_dir) != 'bin':
            asp_dir = join(asp_dir, 'bin')
    else:
        asp_dir = abspath(join('ASP', 'bin'))
    # check ASP install exists
    if not exists(asp_dir):
        raise Exception(f'The following is not the correct path to ASP: {asp_dir}')

    # setup our directory structure
    geoid_dir = join(dirname(in_dir), 'geoid')
    os.makedirs(geoid_dir, exist_ok= True)

    # setup logging
    log_dir = join(geoid_dir, 'logs')
    os.makedirs(log_dir, exist_ok= True)
    log_prefix = 'geoid-tool'
    old_logs = glob(join(log_dir, f'{log_prefix}*.log'))
    if old_logs:
        # gets last run number from the old logs
        vnum = max([int(basename(i).split('.')[0].split('-')[-1].replace('r','')) for i in sorted(old_logs)]) + 1
    else:
        # otherwise sets run number to 1
        vnum = 1

    logging.basicConfig(level=logging.INFO,
                        format=f"(geoid-tool {__name__} %(levelname)s) %(message)s",
                        # saves out to log file and outputs to command line.
                        handlers=[logging.FileHandler(join(log_dir, f'{log_prefix}-r{vnum}.log')),
                                  logging.StreamHandler(sys.stdout)]
                        )
    
    log = logging.getLogger(__name__)

    if debug:
        log.setLevel(logging.DEBUG)  

    # check for white spaces
    if len([i for i in glob(join(in_dir, '*')) if ' ' in i]) > 0:
        log.warning('White spaces found in file paths. Try and remove them?')
        replace_white_spaces(in_dir)

    # Check what kind of file was given by user
    user_file_type = None
    pc_ext = {'.las', '.laz'}
    raster_ext = {'.tif', '.tiff'}
    extension = in_dir.lower().split('.')[-1]
    if f".{extension}" in pc_ext:
        user_file_type = 'Point Cloud'
    elif f".{extension}" in raster_ext:
        user_file_type = 'Raster'
    else:
        raise Exception(f'User data needs to be in following format:{pc_ext} or {raster_ext}')

    # run either transform_pc() or transform_raster()
    if user_file_type == 'Point Cloud':
        log.info(f'Starting run with {user_file_type} data type using PDAL.')
        corrected_data = transform_pc(in_dir=in_dir, 
                                      epsg=epsg, 
                                      asp_dir=asp_dir,
                                      geoid_dir=geoid_dir, 
                                      geoid_model=geoid_model,
                                      log=log)

    elif user_file_type == 'Raster':
        log.info(f'Starting run with {user_file_type} data type using Ames Stereo Pipeline.')
        corrected_data = transform_raster(in_dir=in_dir, 
                                          epsg=epsg, 
                                          asp_dir=asp_dir, 
                                          geoid_dir=geoid_dir,
                                          geoid_model=geoid_model,
                                          transform_command=transform_command,
                                          log=log)

    else:
        raise Exception(f'User data needs to be in following format:{pc_ext} or {raster_ext}')

    # check for success
    if not exists(corrected_data):
        raise Exception('Transformation failed.')
    log.info('Reference data transformed successfully')

    end_time = datetime.now()
    log.info(f"Completed! Run Time: {end_time - start_time}")