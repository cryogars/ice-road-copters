"""
Takes input directory full of .laz (or.las) files and filters+classifies them to DTM laz and DTM tif.

Usage:
    ice-road-pipeline.py <in_dir> [-e user_dem] [-d debug] [-a asp_dir] [-s shp_fp] [-r shp_fp_rfl] [-i imu_data] [-c cal_las] [-k known_rfl] 

Options:
    -e user_dem      Path to user specifed DEM
    -d debug         turns on debugging logging
    -a asp_dir       Directory with ASP binary files
    -s shp_fp        Shapefile to align with
    -g geoid         Is the reference DEM in geoid
    -r shp_fp_rfl    Shapefile to align for reflectance calibration. If given, it is assumed you want grain size output.
                     Additionally, if this mode is selected, the supplied files must be .LAS with extra bytes included with
                     "Intensity as Reflectance" returned by RIEGL.
    -i imu_data      path to helicopter IMU .CSV data used to match data with point cloud using GPS time. 
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
import rasterio
from osgeo import gdal
import numpy as np
import pandas as pd
from datetime import datetime
import logging
import re
import json
import sys
import subprocess
from subprocess import Popen, PIPE
import os

# local imports
from laz2dem import iceroad_logging, las2uncorrectedDEM, cl_call
from laz_align import laz_align
from dir_space_strip import replace_white_spaces
from las_ssa_functions import art_ssa

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

    snow_tif, canopy_tif = laz_align(in_dir = in_dir, align_shp = shp_fp, asp_dir = asp_dir,\
         log = log, input_laz = outlas, canopy_laz = canopy_laz, dem_is_geoid= geoid, las_extra_byte_format=las_extra_byte_format)
    
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

    # mask snow depth from vegetation
    canopyheight = canopyheight.where((canopyheight > snowon + 0.1) | (snowon.isnull()))
    canopyheight.rio.to_raster(canopy_fp)


    ##### START SSA CODE HERE #####
    # due to lidar processing limitations I cannot retain reflectance information using a 
    # transformation. I tested using ASP pc-align and PDAL transformation. Both of these methods
    # I lose the RIEGL reflectance data...
    # SO, instead  deciding that ASP is ran another time (but with translation-only on)
    # And then, we supply X, Y, and Z translation directly to the LAS in lidR package.

    # Command for testing...
    # python ice-road-pipeline.py /Users/brent/Code/ice-road-copters/data/feb9/mcs/ -e /Users/brent/Code/LIDAR/data/QSI_snowfree.tif -a /Users/brent/Code/ice-road-copters/ASP/ -s /Users/brent/Code/ice-road-copters/transform_area/hwy_21/hwy21_utm_edit_v3.shp -r /Users/brent/Code/ice-road-copters/transform_area/Eagle/eagle_res_buffered.shp -i /Users/brent/Code/ice-road-copters/transform_area/Eagle/eagle_res_buffered.shp -d True -c /Users/brent/Code/LIDAR/data/fl_230210_002840/20230209_extraBytes-230210_002840_Scanner_1.las -k 0.142519654

    if shp_fp_rfl:

        # params for SSA function
        k_ice = 1.8983928418833426e-06
        wl = 1064*1e-9
        d_ice = 917
        g = 0.85
        b = 1.6

        # MAKE SURE R IS INSTALLED HERE
        proc = Popen(["which", "R"],stdout=PIPE,stderr=PIPE)
        exit_code = proc.wait()
        if exit_code != 0:
            raise Exception("Please install R on this system and packages: raster, lidR, rlas, dplyr, readr, terra,sf, data.table.")

        # Check whether dem directory exists, if not, make one
        ssa_dir = f'{results_dir}/ssa-calc'
        if not os.path.exists(ssa_dir):
            os.makedirs(ssa_dir)

        # Displaying the parent directory of the script
        scripts_dir = os.path.dirname(__file__)

        # Load in the ASP translational-only matrix from pc-align tool
        asp_matrix_fp = glob(f'{ice_dir}/results/pc-align-translation-only/temp-log*.txt')[0]
        with open(asp_matrix_fp, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if "Translation vector (North-East-Down, meters): Vector3" in line: 
                    n_e_d_shift = re.findall("[+-]?\d+\.\d+", line)
                    n_e_d_shift = ','.join(n_e_d_shift)
                    break

        # get crs
        ref_raster = rasterio.open(snow_tif)
        crs = ref_raster.crs
        pix_size, _ = ref_raster.res
        crs = str(crs).split(":",1)[1]    
        ras_meta = ref_raster.profile

        # Compute slope and aspect from snow-on lidar
        slope_fp = join(ice_dir, f'{basename(in_dir)}-snowon_slope.tif')
        aspect_fp = join(ice_dir, f'{basename(in_dir)}-snowon_aspect.tif')
        os.system(f'gdaldem slope -compute_edges {snow_tif} {slope_fp} -q')
        os.system(f'gdaldem aspect -compute_edges -zero_for_flat {snow_tif} {aspect_fp} -q')
        slope = rio.open_rasterio(slope_fp, masked=True)
        aspect = rio.open_rasterio(aspect_fp, masked=True)

        # Compute surface normal vector from snow-on lidar
        ni_fp = join(ice_dir, f'{basename(in_dir)}-snowon_n_i.tif')
        nj_fp = join(ice_dir, f'{basename(in_dir)}-snowon_n_j.tif')
        nk_fp = join(ice_dir, f'{basename(in_dir)}-snowon_n_k.tif')
        n_i = np.sin(np.radians(aspect)) * np.sin(np.radians(slope))
        n_j = np.cos(np.radians(aspect)) * np.sin(np.radians(slope)) 
        n_k = np.cos(np.radians(slope))
        n_i.rio.to_raster(ni_fp)
        n_j.rio.to_raster(nj_fp)
        n_k.rio.to_raster(nk_fp)        
        
        # Calcs calibration stats from target
        output_csv = f'{ssa_dir}/all-calibration-rfl.csv'
        subprocess.call(["Rscript", 
                          f"{scripts_dir}/las_ssa_cal.r", 
                          cal_las, crs, shp_fp_rfl, n_e_d_shift, output_csv])
        
        # Read in cal data and estimate factor
        df = pd.read_csv(output_csv)
        median_rfl = df['rfl'].median()
        road_cal_factor = known_rfl / median_rfl
        road_cal_factor = str(road_cal_factor) #convert to str for R

        # For each file in <in-dir> 
        for f in os.listdir(in_dir):

            # Rasterize calibrated reflectance and incidence angle
            base_las = os.path.basename(f)
            las_name = os.path.splitext(base_las)[0]
            rfl_fp = f'{ssa_dir}/{las_name}-rfl.las'
            cosi_fp = f'{ssa_dir}/{las_name}-cosi.las'

            # Getting a translated "LAS" file
            # "LAS" in quotations bc I am hiding the cosi and rfl here in "Z"
            # with the intention to do fast IDW in the next step.
            subprocess.call(["Rscript", 
                            f"{scripts_dir}/las_ssa_prep.r", 
                            f, crs, ni_fp, nj_fp, nk_fp, rfl_fp,
                            n_e_d_shift,
                            cosi_fp, road_cal_factor])

            # Run PDAL IDW for rfl_fp and cosi_fp
            # ZACH: resolution is 1.0 in laz2dem.py??
            # Also, it is still fairly slow... not sure if there are any pdal tricks here to speed things up?
            rfl_fp_grid = f'{ssa_dir}/{las_name}-rfl.tif'
            cosi_fp_grid = f'{ssa_dir}/{las_name}-cosi.tif'
            json_path = join(json_dir, 'temp.json')
            json_pipeline = {
                "pipeline": [
                    rfl_fp,
                    {
                        "type":"writers.gdal",
                        "filename":rfl_fp_grid,
                        "resolution":pix_size, 
                        "output_type":"idw"
                    }
                ]
            }

            with open(json_path,'w') as outfile:
                json.dump(json_pipeline, outfile, indent = 2)            
            cl_call(f'pdal pipeline {json_path}', log)      
            os.remove(json_path)

            # Running again for the incidence angle...
            json_pipeline = {
                "pipeline": [
                    cosi_fp,
                    {
                        "type":"writers.gdal",
                        "filename":cosi_fp_grid,
                        "resolution":pix_size, 
                        "output_type":"idw"
                    }
                ]
            }
            with open(json_path,'w') as outfile:
                json.dump(json_pipeline, outfile, indent = 2)
            cl_call(f'pdal pipeline {json_path}', log)        
            os.remove(json_path)
            os.remove(rfl_fp)
            os.remove(cosi_fp)

            # Prepare inputs needed for SSA raster (vectorized operation)
            ssa_fp = f'{ssa_dir}/{las_name}-ssa.tif'
            sini_fp = f'{ssa_dir}/{las_name}-sini.tif'
            theta_fp = f'{ssa_dir}/{las_name}-theta.tif'
            cosi_grid = np.array(gdal.Open(cosi_fp_grid).ReadAsArray())
            ssa_grid = np.ones_like(cosi_grid)
            sini_grid = np.ones_like(cosi_grid)
            sini_grid = np.sin(np.arccos(cosi_grid))
            theta_grid = np.ones_like(cosi_grid)
            theta_grid = np.degrees(np.arccos(-cosi_grid**2 + sini_grid**2 * np.cos(np.radians(180))))
            ref_raster = rasterio.open(cosi_fp_grid)
            ras_meta = ref_raster.profile
            with rasterio.open(ssa_fp, 'w', **ras_meta) as dst:
                 dst.write(ssa_grid, 1)
            with rasterio.open(sini_fp, 'w', **ras_meta) as dst:
                 dst.write(sini_grid, 1)
            with rasterio.open(theta_fp, 'w', **ras_meta) as dst:
                 dst.write(theta_grid, 1)
            ssa_grid = rio.open_rasterio(ssa_fp, masked=True)
            sini_grid = rio.open_rasterio(sini_fp, masked=True)
            cosi_grid = rio.open_rasterio(cosi_fp_grid, masked=True)
            theta_grid = rio.open_rasterio(theta_fp, masked=True)
            rfl_grid = rio.open_rasterio(rfl_fp_grid, masked=True)
            
            # Apply ART - assuming all pixels are snow, directly solve SSA (pretty speedy)
            ssa_grid = (6 * ((4 * np.pi * k_ice) / wl)) / (d_ice * (9*(1-g)) / (16*b) * (-np.log(rfl_grid / ((1.247 + 1.186 * (cosi_grid + cosi_grid) + 5.157 * cosi_grid * cosi_grid + (11.1 * np.exp(-0.087 * theta_grid) + 1.1 * np.exp(-0.014 * theta_grid))) / 4.0 / (cosi_grid + cosi_grid))))**2)

            # Then, go back and clean it up based on CHM and SD
            # For this threshold we use chm is less than 2m
            # .. and snow depth is greater than 8cm
            ssa_grid = ssa_grid.rio.reproject_match(snowdepth)
            ssa_grid = ssa_grid.where((canopyheight <= 2.0) | (snowdepth >= 0.08))
            ssa_grid.rio.to_raster(ssa_fp)     

    end_time = datetime.now()
    log.info(f"Completed! Run Time: {end_time - start_time}")
