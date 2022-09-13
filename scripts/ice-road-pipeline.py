"""
Takes input directory full of .laz files and filters+classifies them to DTM laz and DTM tif.

Usage:
    ice-road-pipeline.py <in_dir> [-e user_dem] [-d debug] [-a asp_dir] [-s shp_fp]

Options:
    -e user_dem      Path to user specifed DEM
    -d debug         turns on debugging logging
    -a asp_dir       Directory with ASP binary files
    -s shp_fp        Shapefile to align with
"""

from docopt import docopt
from glob import glob
from os.path import abspath, join, basename, isdir
from laz_align import laz_align
import rioxarray as rio
from datetime import datetime
import logging
import sys
import os

# local imports
from laz2dem import iceroad_logging, las2uncorrectedDEM
from laz_align import laz_align
from dir_space_strip import replace_white_spaces


if __name__ == '__main__':
    start_time = datetime.now()
    # get command line args
    args = docopt(__doc__)
    user_dem = args.get('-e')
    if user_dem:
        user_dem = abspath(user_dem)
    debug = args.get('-d')
    asp_dir = args.get('-a')
    if asp_dir:
        asp_dir = abspath(asp_dir)
    else:
        asp_dir = abspath(join('ASP', 'bin'))
    shp_fp = args.get('-s')
    if shp_fp:
        shp_fp = abspath(shp_fp)
        if isdir(shp_fp):
            raise Exception("Provide to .shp file to use. Not directory.")
    else:
        raise Exception("Provide filepath to .shp file for alignment with -s flag.")

    in_dir = args.get('<in_dir>')
    # convert to abspath
    in_dir = abspath(in_dir)
    # setup our directory structure
    ice_dir = join(in_dir, 'ice-road')
    os.makedirs(ice_dir, exist_ok= True)
    results_dir = join(ice_dir, 'results')
    os.makedirs(results_dir, exist_ok= True)
    json_dir = join(ice_dir, 'json')
    os.makedirs(results_dir, exist_ok= True)

    # setup logging
    log_dir = join(ice_dir, 'logs')
    os.makedirs(log_dir, exist_ok= True)
    log_prefix = 'ice-road-pipeline'
    old_logs = glob(join(log_dir, f'{log_prefix}*.log'))
    if old_logs:
        # gets last run number from the old logs
        vnum = max([int(basename(i).split('.')[0].split('-')[-1].replace('r','')) for i in sorted(old_logs)]) + 1
    else:
        # otherwise sets run number to 1
        vnum = 1

    logging.basicConfig(level=logging.INFO,
    format=f"(ice-road-copters {__name__} %(levelname)s) %(message)s",
    # saves out to log file and outputs to command line.
    handlers=[
        logging.FileHandler(join(log_dir, f'{log_prefix}-r{vnum}.log')),
        logging.StreamHandler(sys.stdout)]
    )
    log = logging.getLogger(__name__)
    if debug:
        log.setLevel(logging.DEBUG)  

    # run main functions
    log.info('Starting laz2uncorrectedDEM')
    log.info(f'Using in_dir: {in_dir}, user_dem: {user_dem}')
    outtif, outlas, canopy_laz = las2uncorrectedDEM(in_dir, debug, log, user_dem = user_dem)
    
    log.info('Starting ASP laz align')
    log.info(f'Using in_dir: {in_dir}, shapefile: {shp_fp}, ASP dir: {asp_dir}')

    aligned_tif = laz_align(work_dir = in_dir, align_shp = shp_fp, asp_dir = asp_dir,\
         log = log, input_laz = outlas, canopy_laz = canopy_laz)
    if aligned_tif == -1:
        raise Exception('Failed to align to shapefile.')
    
    # difference two rasters to find snow depth
    ref_dem_path = join(in_dir, 'results/ref_PC.tif')
    snow_dem_path = join(in_dir, 'results/pc-grid/run-DEM.tif')
    snow_depth_path = join(in_dir, 'results/snowdepth.tif')
    snowoff = rio.open_rasterio(ref_dem_path, masked=True)
    snowon = rio.open_rasterio(snow_dem_path, masked=True) 
    snowon_matched = snowon.rio.reproject_match(snowoff)
    snowdepth = snowon_matched - snowoff
    snowdepth.rio.to_raster(snow_depth_path)

    end_time = datetime.now()
    log.info(f"Completed! Run Time: {end_time - start_time}")
