"""
Takes input directory full of .laz (or.las) files and filters+classifies them to DTM laz and DTM tif.

Usage:
    ice-road-pipeline.py <in_dir> [-e user_dem] [-d debug] [-a asp_dir] [-s shp_fp] [-b buffer_meters] [-r shp_fp_rfl] [-i imu_data] [-c cal_las] [-k known_rfl] 

Options:
    -e user_dem      Path to user specifed DEM
    -d debug         turns on debugging logging
    -a asp_dir       Directory with ASP binary files
    -s shp_fp        Shapefile to align with
    -b buffer_meters Total width for the transform area
    -g geoid         Is the reference DEM in geoid
    -r shp_fp_rfl    Shapefile to align for reflectance calibration. If given, it is assumed you want grain size output.
                     Additionally, if this mode is selected, the supplied files must be .LAS with extra bytes included with
                     "Intensity as Reflectance" returned by RIEGL.
    -i imu_data      path to helicopter IMU .CSV or.TXT data used to match data with point cloud using GPS time.
                     Column names must include ['Time[s]', 'Easting[m]', 'Northing[m]', 'Height[m]'] 
    -c cal_las       path to .LAS used for calibration of the apparent reflectance for 1064nm of lidar sensor.
                     To avoid confusion, please supply this file in a different directory from <in_dir>.
    -k known_rfl     Known intrinsic reflectance at 1064nm (float/real) for target identified in shp_fp_rfl.

"""

from cmath import exp
from docopt import docopt
from glob import glob
from os.path import abspath, join, basename, isdir
from laz_align import laz_align
import rioxarray as rio
from rasterio.crs import CRS
from datetime import datetime
import logging
import sys
import os

# local imports
from laz2dem import iceroad_logging, las2uncorrectedDEM, cl_call
from laz_align import laz_align
from dir_space_strip import replace_white_spaces
from las2grain import ssa_pipeline

if __name__ == '__main__':
    start_time = datetime.now()
    # get command line args
    args = docopt(__doc__)
    user_dem = args.get('-e')
    if user_dem:
        user_dem = abspath(user_dem)
    geoid = args.get('-g')
    if not user_dem:
        geoid = True
    debug = args.get('-d')
    asp_dir = args.get('-a')
    if asp_dir:
        asp_dir = abspath(asp_dir)
        if basename(asp_dir) != 'bin':
            asp_dir = join(asp_dir, 'bin')
    else:
        asp_dir = abspath(join('ASP', 'bin'))
    shp_fp = args.get('-s')
    if shp_fp:
        shp_fp = abspath(shp_fp)
        if isdir(shp_fp):
            raise Exception("Provide to .shp file to use. Not directory.")
        elif not shp_fp.endswith('.shp'):
            raise Exception("Provide fp to .shp file to use.")
    else:
        raise Exception("Provide filepath to .shp file for alignment with -s flag.")
    
    buffer_meters = args.get('-b')
    if buffer_meters:
        known_rfl = float(buffer_meters)
    else:
        buffer_meters = 3.0

    shp_fp_rfl = args.get('-r')
    if shp_fp_rfl:
        shp_fp_rfl = abspath(shp_fp_rfl)
        las_extra_byte_format = True
    else:
        las_extra_byte_format = False

    imu_data = args.get('-i')
    if imu_data:
        imu_data = abspath(imu_data)

    cal_las = args.get('-c')
    if cal_las:
        cal_las = abspath(cal_las)

    known_rfl = args.get('-k')
    if known_rfl:
        known_rfl = float(known_rfl)

    in_dir = args.get('<in_dir>')
    # convert to abspath
    in_dir = abspath(in_dir)
    
    
    # setup our directory structure
    ice_dir = join(in_dir, 'ice-road')
    os.makedirs(ice_dir, exist_ok= True)
    results_dir = join(ice_dir, 'results')
    os.makedirs(results_dir, exist_ok= True)
    json_dir = join(ice_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)

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

    # check for white spaces
    if len([i for i in glob(join(in_dir, '*')) if ' ' in i]) > 0:
        log.warning('White spaces found in file paths. Try and remove them?')
        replace_white_spaces(in_dir)
        # raise Exception('File paths contains spaces. Please remove with the dir_space_strip.py script.')

    # run main functions
    log.info('Starting laz2uncorrectedDEM')
    log.info(f'Using in_dir: {in_dir}, user_dem: {user_dem}')
    outtif, outlas, canopy_laz = las2uncorrectedDEM(in_dir, debug, log, 
                                                    user_dem = user_dem, 
                                                    las_extra_byte_format = las_extra_byte_format)
    
    log.info('Starting ASP laz align')
    log.info(f'Using in_dir: {in_dir}, shapefile: {shp_fp}, ASP dir: {asp_dir}')

    snow_tif, canopy_tif = laz_align(in_dir = in_dir, align_shp = shp_fp, 
                                     asp_dir = asp_dir,log = log, input_laz = outlas, 
                                     canopy_laz = canopy_laz, dem_is_geoid= geoid, 
                                     buffer_meters=buffer_meters,
                                     las_extra_byte_format=las_extra_byte_format)
    
    # clean up after ASP a bit
    for fp in os.listdir(ice_dir):
        if fp.endswith(".txt"):
            os.remove(join(ice_dir, fp))
        if fp.endswith('-DEM.tif'):
            os.rename(join(ice_dir, fp), join(ice_dir, fp.replace('-DEM','')))
    snow_tif = snow_tif.replace('-DEM','')
    canopy_tif = canopy_tif.replace('-DEM','')
    
    # difference two rasters to find snow depth
    ref_dem_path = join(results_dir, 'dem.tif')
    snow_depth_path = join(ice_dir, f'{basename(in_dir)}-snowdepth.tif')
    snowoff = rio.open_rasterio(ref_dem_path, masked=True)
    snowon = rio.open_rasterio(snow_tif, masked=True) 
    snowon_matched = snowon.rio.reproject_match(snowoff)
    snowdepth = snowon_matched - snowoff

    # Write reproject-match SnowDepth and SnowOn to disk
    snowdepth.rio.to_raster(snow_depth_path)
    snowon_matched.rio.to_raster(snow_tif)

    # difference two rasters to find canopy height
    ref_dem_path = join(results_dir, 'dem.tif')
    canopy_fp = join(ice_dir, f'{basename(in_dir)}-canopyheight.tif')
    canopy = rio.open_rasterio(canopy_tif, masked=True) 
    matched = canopy.rio.reproject_match(snowoff)
    canopyheight = matched - snowoff

    # mask snow depth from vegetation (null=0)
    canopyheight = canopyheight.where((canopyheight > snowdepth + 0.1) | (snowdepth.isnull()), other=0)
    canopyheight.rio.to_raster(canopy_fp)

    # estimate ssa from reflectance using AART at 1064 nm.
    if shp_fp_rfl:
        ssa_grid = ssa_pipeline(snowon_matched, cal_las, 
                                shp_fp_rfl, imu_data, known_rfl, 
                                json_dir, results_dir, ice_dir, in_dir, snow_tif)
        ssa_fp = join(ice_dir, f'{basename(in_dir)}-ssa.tif')

        # Set threshold for snow depth and canopy height
        ssamask = canopyheight.rio.reproject_match(ssa_grid)
        snowdepth = snowdepth.rio.reproject_match(ssa_grid)
        ssamask = ssamask.where(ssamask <= 3, -100)
        ssa_grid = ssa_grid.where(ssamask != -100)
        ssa_grid = ssa_grid.where(snowdepth >= 0.08)
        ssa_grid.rio.to_raster(ssa_fp)

    end_time = datetime.now()
    log.info(f"Completed! Run Time: {end_time - start_time}")