"""
Takes input directory full of .laz files and filters+classifies them to DTM laz and DTM tif.

Usage:
    laz2dem.py <in_dir> [-d debug]

Options:
    -d debug      turns on debugging logging  [default: True]
"""
import json
import logging
import os
import shlex
import subprocess
import sys
from datetime import datetime
from glob import glob
from os.path import abspath, basename, dirname, exists, isdir, join, expanduser

import laspy
import py3dep
import pyproj
import rioxarray as rxa
from docopt import docopt
from rasterio.enums import Resampling
from shapely.geometry import box
from shapely.ops import transform

log = logging.getLogger(__name__)

def cl_call(command, log):
    """
    Runs shell commands in python and returns output
    Got this from a stack overflow but can't find it now...

    Parameters:
    command (str or list): list of commands. if string is passed we will try and
    parse to list using shelex
    """
    if type(command) == str:
        command = shlex.split(command)
    log.info('Subprocess: "' + ' '.join(command) + '"')

    process = subprocess.Popen(command, 
                           stdout=subprocess.PIPE,
                           universal_newlines=True)

    while True:
        output = process.stdout.readline()
        log.info(output.strip())
        # Do something else
        return_code = process.poll()
        if return_code is not None:
            log.info(f'RETURN CODE {return_code}')
            # Process has finished, read rest of the output 
            for output in process.stdout.readlines():
                log.info(output.strip())
            break

def create_json_pipeline(in_fp, outlas, outtif, dem_fp, json_name = 'las2dem', json_dir = './jsons'):
    """
    Creates JSON Pipeline for standard las point cloud to DTM.
    Filters include: dem, elm, outlier
    SMRF Classifier and writes ground classified points to las and tif

    Parameters:
    in_fp (str): filepath to las file to be run
    outlas (str): filepath to save dtm las
    outtif (str): filepath to save dtm tif
    json_name (str) [optional]: name of json to save [default: las2dem.json]
    json_dir (str) [optional]: name of json subdirectory to create [default: ./json]

    Returns:
    json_to_use (str): filepath of created json pipeline
    """
    ## make sure path is in format pdal likes
    in_fp = abspath(in_fp)
    outlas = abspath(outlas)
    outtif = abspath(outtif)

    assert exists(in_fp), f'In filepath {in_fp} does not exist'
    assert exists(dem_fp), f'DEM filepath {in_fp} does not exist'

    # good docs on types of filters used: https://pdal.io/stages/filters.html#ground-unclassified
    # Reads in mosaiced las file
    reader = {"type": "readers.las", "filename": in_fp
    }

    # Filters out points with 0 returns
    mongo_filter = {"type": "filters.mongo",\
        "expression": {"$and": [\
            {"ReturnNumber": {"$gt": 0}},\
                {"NumberOfReturns": {"$gt": 0}} ] } 
    }
    # Filter out points far away from our dem
    dem_filter = {
            "type":"filters.dem",
            "raster":dem_fp,
            "limits":"Z[25:35]"
    }
    # Extended Local Minimum filter
    elm_filter = {"type": "filters.elm"
    }
    # Outlier filter
    outlier_filter = {"type": "filters.outlier",\
        "method": "statistical",\
            "mean_k": 12,\
                "multiplier": 2.2
    }
    # SMRF classifier for ground
    smrf_classifier = {"type": "filters.smrf",\
        "ignore": "Classification[7:7], NumberOfReturns[0:0], ReturnNumber[0:0]"
    }
    # Select ground points only
    smrf_selecter = { 
            "type":"filters.range",
            "limits":"Classification[2:2]"
    }
    # Write las file
    las_writer = {"type": "writers.las",\
#     "where": "Classification[2:2]",\
        "filename":outlas
    }
    # Write tif file
    tif_writer = {"type": "writers.gdal",\
    #     "where": "Classification[2:2]",\
            "filename":outtif,
            "resolution":1.0,
            "output_type":"idw"
    }
    # set up pipeline
    pipeline = [reader, mongo_filter, dem_filter, elm_filter, outlier_filter, smrf_classifier,smrf_selecter, las_writer, tif_writer]
    # make json dir and fp
    log.debug(f"Making JSON dir at {json_dir}")
    os.makedirs(json_dir, exist_ok= True)
    json_name = json_name.replace('.json','')
    json_to_use = join(json_dir, f'{json_name}.json')
    # write json fp out
    with open(json_to_use,'w') as outfile:
        json.dump(pipeline, outfile, indent = 2)
    # add logging message for success #

    return json_to_use

def mosaic_laz(in_dir, log, out_fp = 'merge.laz', laz_prefix = ''):
    """
    Generates and run PDAL mosaic command.

    Parameters:
    in_dir (str): fp to directory full of .laz files to mosaic
    out_fp (str) [optional]: out filepath to save [default: ./merge.laz]
    laz_prefix (str) [optional]: prefix to append in case there are .laz files 
    to avoid mosaicing [default: ""]
    Returns:
    mosaic_fp (str): filepath to mosaic output file
    """
    assert isdir(in_dir), f'{in_dir} is not a directory'
    # generate searching command
    in_str = ' '.join(glob(join(in_dir, f'{laz_prefix}*.laz')))
    # out fp to save to
    mosaic_fp = join(in_dir, out_fp)
    # set up mosaic command
    mosaic_cmd = f'pdal merge {in_str} {mosaic_fp}'
    log.debug(f"Using mosaic command: {mosaic_cmd}")
    # run mosaic command
    cl_call(mosaic_cmd, log)
    
    return mosaic_fp

def download_dem(las_fp, dem_fp = 'dem.tif'):
    """
    Reads the crs and bounds of a las file and downloads a DEM from py3dep
    Must be in the CONUS.

    Parameters:
    las_fp (str): filepath to las file to get bounds and crs
    dem_fp (str) [optional]: filepath to save DEM at. [default = './dem.tif']

    Returns:
    crs (pyproj CRS): CRS object from las header
    project (shapely transform): shapely transform used in conversion
    """
    # read crs of las file
    with laspy.open(las_fp) as las:
        hdr = las.header
        crs = hdr.parse_crs()
    log.debug(f"CRS used is {crs}")
    # create transform from wgs84 to las crs
    wgs84 = pyproj.CRS('EPSG:4326')
    project = pyproj.Transformer.from_crs(crs, wgs84 , always_xy=True).transform
    # calculate bounds of las file in wgs84
    utm_bounds = box(hdr.mins[0], hdr.mins[1], hdr.maxs[0], hdr.maxs[1])
    wgs84_bounds = transform(project, utm_bounds)
    # download dem inside bounds
    dem_wgs = py3dep.get_map('DEM', wgs84_bounds, resolution=1, crs='EPSG:4326')
    log.debug(f"DEM bounds: {dem_wgs.rio.bounds()}. Size: {dem_wgs.size}")
    # reproject to las crs and save
    dem_utm = dem_wgs.rio.reproject(crs, resampling = Resampling.cubic_spline)
    dem_utm.rio.to_raster(dem_fp)
    log.debug(f"Saved to {dem_fp}")
    return dem_fp, crs, project

def las2uncorrectedDEM(in_dir, debug, log, user_dem = None):
    """
    Takes a input directory of laz files. Mosaics them, downloads DEM within their bounds,
    builds JSON pipeline, and runs PDAL pipeline of filter, classifying and saving DTM.

    Parameters:
    in_dir (str): filepath to directory to run in
    debug (bool): lots of yakety yak or not?

    Returns:
    outtif (str): filepath to output DTM tiff
    outlas (str): filepath to output DTM laz file
    """
    # log_dir = join(in_dir, 'logs')
    # log = iceroad_logging(log_dir, debug, log_prefix = 'filter_classify')
    #set start time
    # start_time = datetime.now()
    # checks on directory and user update
    assert isdir(in_dir), f'Provided: {in_dir} is not a directory. Provide directory with .laz files.'
    log.info(f"Working in directory: {in_dir}")
    os.chdir(in_dir)

    # set up sub directories
    results_dir = join(in_dir, 'results')
    os.makedirs(results_dir, exist_ok= True)

    # check for overwrite
    outtif = join(results_dir, f'{basename(in_dir)}.tif')
    outlas = join(results_dir, f'{basename(in_dir)}.laz')
    if exists(outtif):
        while True:
            ans = input("Uncorrected tif already exists. Enter y to overwrite and n to use existing:")
            if ans.lower() == 'n':
                return outtif, outlas
            elif ans.lower() == 'y':
                break

    # mosaic
    log.info("Starting to mosaic las files...")
    las_fps = glob(join(in_dir, '*.laz'))
    log.debug(f"Number of las files: {len(las_fps)}")
    mosaic_fp = join(results_dir, 'merge.laz')
    mosaic_fp = mosaic_laz(in_dir, out_fp=mosaic_fp, log = log)

    if not exists(mosaic_fp):
        log.warning('No mosaic created')
        return -1

    # Allowing the code to use user input DEM
    dem_fp = join(results_dir, 'dem.tif')
    if user_dem is False:
        log.info("Starting DEM download...")
        _, crs, project = download_dem(mosaic_fp, dem_fp = dem_fp )
        log.debug(f"Downloaded dem to {dem_fp}")
    else:
        log.info("User DEM specified. Skipping DEM download...")
        cl_call('cp '+ user_dem +' '+ dem_fp, log) #to ensure const. filenames for next step
    if not exists(join(results_dir, 'dem.tif')):
        log.warning('No DEM downloaded')
        return -1

    log.info("Creating JSON Pipeline...")
    json_to_use = create_json_pipeline(in_fp = mosaic_fp, outlas = outlas, outtif = outtif, dem_fp = dem_fp)
    log.debug(f"JSON to use is {json_to_use}")

    log.info("Running PDAL pipeline")
    if debug:
        pipeline_cmd = f'pdal pipeline -i {json_to_use} -v 8'
    else:
        pipeline_cmd = f'pdal pipeline -i {json_to_use}'
    cl_call(pipeline_cmd, log)

    # end_time = datetime.now()
    # log.info(f"Completed! Run Time: {end_time - start_time}")

    return outtif, outlas

def iceroad_logging(log_dir, debug, log_prefix = 'las2uncorrectedDEM' ):
    os.makedirs(log_dir, exist_ok= True)

    old_logs = glob(join(log_dir, f'{log_prefix}*.log'))
    if old_logs:
        vnum = max([int(basename(i).split('.')[0].split('-')[-1].replace('r','')) for i in sorted(old_logs)]) + 1
    else:
        vnum = 1

    logging.basicConfig(level=logging.INFO,
    format=f"(ice-road-copters {__name__} %(levelname)s) %(message)s",
    handlers=[
        logging.FileHandler(join(log_dir, f'{log_prefix}-r{vnum}.log')),
        logging.StreamHandler(sys.stdout)]
    )
    log = logging.getLogger(__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    
    return log

if __name__ == '__main__':
    # get command line args
    args = docopt(__doc__)
    debug = args.get('-d')
    in_dir = args.get('<in_dir>')
    # convert to abspath
    in_dir = abspath(in_dir)
    # setup logging

    # run main function
    outtif, outlas = las2uncorrectedDEM(in_dir, debug, log)
