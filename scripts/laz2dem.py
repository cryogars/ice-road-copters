"""
Takes input directory full of .laz files and filters+classifies them to DTM laz and DTM tif.

Usage:
    laz2dem.py <in_dir> [-d debug] [-S smrf_scalar] [-L smrf_slope] [-T smrf_threshold] [-W smrf_window]

Options:
    -d debug          turns on debugging logging  [default: True]
    -S smrf_scalar    (Optional) Override SMRF scalar parameter (float)
    -L smrf_slope     (Optional) Override SMRF slope parameter (float)
    -T smrf_threshold (Optional) Override SMRF threshold parameter (float)
    -W smrf_window    (Optional) Override SMRF window size (float)
"""
import json
import logging
import os
import shlex
import subprocess
import sys
from datetime import datetime
from typing import Optional
from glob import glob
from os.path import abspath, basename, exists, isdir, join

import laspy
import py3dep
import pyproj
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

def create_json_pipeline(
        in_fp: str, outlas: str, outtif: str,
        json_name: Optional[str] = 'las2unaligned',
        json_dir: Optional[str] = './jsons',
        canopy: bool = False,
        smrf_overrides: dict[str, float] | None = None
) -> str:
    """
    Creates JSON Pipeline for standard las point cloud to DTM.
    Filters include: elm, outlier
    SMRF Classifier and writes ground classified points to las and tif

    Parameters:
    in_fp: filepath to las file to be run
    outlas: filepath to save dtm las
    outtif: filepath to save dtm tif
    json_name: name of json to save [default: las2dem.json]
    json_dir: name of json subdirectory to create [default: ./json]
    smrf_overrides: optional SMRF parameter (scaler, slope threshold, and window) overrides.

    Returns:
    json_to_use: filepath of created json pipeline
    """
    ## make sure path is in format pdal likes
    in_fp = abspath(in_fp)
    outlas = abspath(outlas)
    outtif = abspath(outtif)

    assert exists(in_fp), f'In filepath {in_fp} does not exist'

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
    smrf_classifier = {"type": "filters.smrf",
        "ignore": "Classification[7:7], NumberOfReturns[0:0], ReturnNumber[0:0]"
    }

    if smrf_overrides:
        smrf_classifier.update(smrf_overrides)
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

    first_returns = {"type": "filters.range",\
            "limits":"returnnumber[1:1]"
    }

    # set up pipeline
    if canopy:
        pipeline = [reader, first_returns, las_writer]
    else:
        pipeline = [
            reader,
            mongo_filter,
            elm_filter,
            outlier_filter,
            smrf_classifier,
            smrf_selecter,
            las_writer,
            tif_writer
        ]
        
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

def mosaic_laz(in_dir, las_extra_byte_format, log, out_fp = 'unaligned_merged.laz', laz_prefix = ''):
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
    if las_extra_byte_format is True:
        in_str = ' '.join(glob(join(in_dir, f'{laz_prefix}*.las')))
    else:
        in_str = ' '.join(glob(join(in_dir, f'{laz_prefix}*.laz')))
    # out fp to save to
    mosaic_fp = join(in_dir, out_fp)
    # set up mosaic command
    mosaic_cmd = f'pdal merge {in_str} {mosaic_fp}'
    log.debug(f"Using mosaic command: {mosaic_cmd}")
    # run mosaic command
    cl_call(mosaic_cmd, log)
    
    return mosaic_fp

def las2uncorrectedDEM(
        in_dir: str,
        log: logging.Logger,
        debug: bool = False,
        las_extra_byte_format: bool = False,
        smrf_overrides: dict[str, float] | None = None,
        force_overwrite: bool = False
) -> tuple[str, str, str]:
    """
    Takes a input directory of laz files. Mosaics them, builds JSON pipeline,
    and runs PDAL pipeline of filter, classifying and saving DTM.

    Parameters:
    -----------
    in_dir: filepath to directory to run in
    log: logger instance
    debug (default is False): lots of yakety yak or not? 
    smrf_overrides: optional SMRF parameter (scaler, slope threshold, and window) overrides.

    Returns:
    --------
    Filepaths to generated DTM TIF, DTM LAS/LAZ, and canopy LAS.
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
    ice_dir = join(in_dir, 'ice-road')
    os.makedirs(ice_dir, exist_ok= True)
    results_dir = join(ice_dir, 'results')
    os.makedirs(results_dir, exist_ok= True)
    json_dir =  join(ice_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)

    # check for overwrite
    outtif = join(results_dir, f'{basename(in_dir)}_unaligned.tif')
    outlas = join(results_dir, f'{basename(in_dir)}_unaligned.laz')
    canopy_laz = join(results_dir, f'{basename(in_dir)}_canopy_unaligned.laz')
    
    if exists(outtif) and not force_overwrite:
        log.info("Reusing existing uncorrected TIF")
        return outtif, outlas, canopy_laz
    elif exists(outtif) and force_overwrite:
        log.info("Overwriting existing uncorrected TIF")

    # mosaic
    log.info("Starting to mosaic las files...")
    if las_extra_byte_format is True:
        las_fps = glob(join(in_dir, '*.las'))
    else:
        las_fps = glob(join(in_dir, '*.laz'))
    log.debug(f"Number of las files: {len(las_fps)}")
    mosaic_fp = join(results_dir, 'unfiltered_merge.laz')
    mosaic_fp = mosaic_laz(in_dir,las_extra_byte_format, out_fp=mosaic_fp, log = log)

    if not exists(mosaic_fp):
        log.warning('No mosaic created')
        return -1

    # DTM creation
    log.info("Creating DTM Pipeline...")
    json_to_use = create_json_pipeline(
        in_fp=mosaic_fp,
        outlas=outlas,
        outtif=outtif,
        json_dir=json_dir,
        smrf_overrides=smrf_overrides
    )
    log.debug(f"JSON to use is {json_to_use}")

    log.info("Running DTM pipeline")
    if debug:
        pipeline_cmd = f'pdal pipeline -i {json_to_use} -v 8'
    else:
        pipeline_cmd = f'pdal pipeline -i {json_to_use}'
    cl_call(pipeline_cmd, log)

    # DSM creation
    log.info("Creating Canopy Pipeline...")
    json_to_use = create_json_pipeline(
        in_fp=mosaic_fp,
        outlas=canopy_laz,
        outtif=canopy_laz.replace('laz', 'tif'),
        json_dir=json_dir,
        canopy=True,
        json_name='canopy',
        smrf_overrides=smrf_overrides
    )
    log.debug(f"JSON to use is {json_to_use}")

    log.info("Running Canopy pipeline")
    if debug:
        pipeline_cmd = f'pdal pipeline -i {json_to_use} -v 8'
    else:
        pipeline_cmd = f'pdal pipeline -i {json_to_use}'
    cl_call(pipeline_cmd, log)


    # end_time = datetime.now()
    # log.info(f"Completed! Run Time: {end_time - start_time}")

    return outtif, outlas, canopy_laz


def filter_dem(
    input_laz: str,
    dem_fp: str,
    outlas: str,
    json_dir: str,
    log: logging.Logger
) -> str:
    """
    Apply PDAL DEM filtering to a point cloud (https://pdal.io/en/stable/stages/filters.dem.html).

    Parameters
    ----------
    input_laz: path to a LAS/LAZ file to filter
    dem_fp: filepath to a DEM raster used for filtering
    outlas: output path for the DEM-filtered LAS/LAZ file
    json_dir: directory path to store the generated PDAL JSON pipeline
    log: logger instance

    Returns
    -------
    File path to the generated DEM-filtered LAS/LAZ.
    """

    input_laz = abspath(input_laz)
    outlas = abspath(outlas)
    dem_fp = abspath(dem_fp)

    assert exists(input_laz), f"Input LAS/LAZ not found: {input_laz}"
    assert exists(dem_fp), f"DEM raster not found: {dem_fp}"

    log.debug(f"Using {dem_fp} raster for filtering")

    pipeline = [
        {"type": "readers.las", "filename": input_laz},
        {"type": "filters.dem", "raster": dem_fp, "limits": "Z[-10:10]"},
        {"type": "writers.las", "filename": outlas}
    ]

    os.makedirs(json_dir, exist_ok=True)
    json_to_use = join(json_dir, "filter_dem.json")

    with open(json_to_use, "w") as outfile:
        json.dump(pipeline, outfile, indent=2)

    log.info(f"Running DEM filter on {input_laz}")
    cl_call(f"pdal pipeline {json_to_use}", log)

    if not exists(outlas):
        raise RuntimeError(f"filter_dem(): failed to produce output file: {outlas}")
    
    log.info(f"DEM filtering complete: {outlas}")
    return outlas


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
    
    smrf_overrides = {}
    for flag, key in [('-S', 'scalar'), ('-L', 'slope'), ('-T', 'threshold'), ('-W', 'window')]:
        value = args.get(flag)
        if value is not None:
            smrf_overrides[key] = float(value)

    # run main function
    outtif, outlas, _ = las2uncorrectedDEM(
        in_dir=in_dir,
        debug=debug,
        log=log,
        user_dem=None,
        las_extra_byte_format=False,
        smrf_overrides=smrf_overrides
    )
