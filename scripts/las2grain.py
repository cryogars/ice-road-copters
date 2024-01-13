from glob import glob
from os.path import join, basename
import rioxarray as rio
import rasterio
import numpy as np
import pandas as pd
import logging
import re
import json
import subprocess
from subprocess import Popen, PIPE
import os

# local imports
from laz2dem import cl_call


log = logging.getLogger(__name__)


def normalized_reflectance(snow_tif, cal_las, shp_fp_rfl,
                           imu_data, known_rfl, json_dir,
                           results_dir, ice_dir, in_dir):
    '''
    This function utilizes lidR package to leverage the use of key features such as translation,
    sample grid data, and apply math operations. With this being an R package, this function requires
    the installation of R and the following packages:

    raster, lidR, rlas, dplyr, readr, terra,sf, data.table

    In summary,
    1 of the n LAS fight data is input, LAS is shifted in X, Y, and Z direction based on ASP pc_align with translation only,
    Heli IMU data (X,Y,Z,t) is used to build vector between surface and sensor --> used to estimate incidence angles, 
    Reflectance is normalized + cleaned, finally, a new .LAS file is created where X and Y are position, but Z is the reflectance.

    Currently it is called within the SSA function below `ssa_pipeline`.

    '''

    # Make sure R is installed
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
    # NOTE: this is assuming translation file doesn't change for a dataset.
    asp_matrix_fp = glob(f'{ice_dir}/results/pc-align-translation-only/temp-log*.txt')[0]
    with open(asp_matrix_fp, 'r') as f:
        lines = f.readlines()
        for line in lines:
            asp_line_found=0
            if "Translation vector (North-East-Down, meters): Vector3" in line: 
                n_e_d_shift = re.findall("[+-]?\d+\.\d+", line)
                n_e_d_shift = ','.join(n_e_d_shift)
                asp_line_found=1
                log.info(f'North, East, Down shift: {n_e_d_shift}.')
                break
    if asp_line_found==0:
        raise Exception('There was an issue finding translation vector in the text file. ASP pc-align tool may have changed their output if this is the case.')
    
    # get crs
    ref_raster = rasterio.open(snow_tif)
    crs = ref_raster.crs
    pix_size, _ = ref_raster.res
    crs = str(crs).split(":",1)[1]
    log.info(f'Using pixel size: {pix_size} m and using CRS: {crs}.') 

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
    log.info(f'Surface normal computed.')       
    
    # Calcs calibration stats from target
    output_csv = f'{ssa_dir}/all-calibration-rfl.csv'
    subprocess.call(["Rscript", f"{scripts_dir}/las_ssa_cal.r", 
                     cal_las, crs, shp_fp_rfl, n_e_d_shift, 
                     output_csv, imu_data, str(pix_size)], 
                     stdout=log, stderr=log)
    
    # Read in cal data and estimate factor
    df = pd.read_csv(output_csv)
    median_rfl = df['rfl'].median()
    road_cal_factor = known_rfl / median_rfl
    road_cal_factor = str(road_cal_factor) #convert to str for R

    # For each file in <in-dir> 
    for f in glob(in_dir+'/*.las'):   

        # Set names
        base_las = os.path.basename(f)
        las_name = os.path.splitext(base_las)[0]

        # Quick check to make sure not using wrong files
        if las_name == ".DS_Store" or las_name == "ice-road":
            continue 
        
        # Set file path to this specific .las in loop
        rfl_fp = f'{ssa_dir}/{las_name}-rfl.las'

        # Getting a translated "LAS" file
        # "LAS" in quotations bc I am hiding the rfl here in "Z"
        # with the intention to do fast IDW in the next step.
        subprocess.call(["Rscript", 
                         f"{scripts_dir}/las_ssa_prep.r", 
                         f, crs, ni_fp, nj_fp, nk_fp, rfl_fp,
                         n_e_d_shift,road_cal_factor,imu_data], 
                         stdout=log, stderr=log)

    rfl_fp = f'{ssa_dir}/rfl-merged.las'
    rfl_fp_grid = f'{ssa_dir}/rfl-merged.tif'

    # Prep PDAL commands
    json_path = join(json_dir, 'combine-rfl.json')
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

    # Run PDAL commands
    inp_str = ' '.join(glob(join(ssa_dir, '*.las')))
    cl_call(f'pdal merge {inp_str} {rfl_fp}', log)       
    cl_call(f'pdal pipeline {json_path}', log)  

    return rfl_fp_grid




def aart_1064(rfl_grid, cosi=1, g=0.85, b=1.6):

    '''
    Solves SSA for AART at 1064 nm wavelength. 
    Assumes pixel is 100% snow and that lidar is pure backscatter.
    g=0.85 and b=1.6 are from Libois et al. (2013), however, other shapes may be used.

    The cosine of light incidence angle with surface normal (cosi) 
    is normalized by incidence angle in las_ssa_prep.R, cosi is fixed at 1.

    '''

    log.info(f'Computing SSA with AART using g={g} and b={b} .')  

    theta = 180 # based on data it is almost always 179-180

    # This is the imaginary refractive index of ice for the wavelength at 1064 nm (fixed). 
    k_ice = 1.8983928418833426e-06

    # This is the wavelength in nanometers, then multiplied to meters for computation (fixed).
    wl = 1064*1e-9

    # This is the density of ice in kg/m3 (fixed).
    d_ice = 917

    # Run AART - assuming all pixels are snow, directly solve SSA
    ssa_grid = (6 * ((4 * np.pi * k_ice) / wl)) / (d_ice * (9*(1-g)) / (16*b) * (-np.log(rfl_grid / ((1.247 + 1.186 * (cosi + cosi) + 5.157 * cosi * cosi + (11.1 * np.exp(-0.087 * theta) + 1.1 * np.exp(-0.014 * theta))) / 4.0 / (cosi + cosi))))**2)

    return ssa_grid





def ssa_pipeline(snow_tif, snowoff, cal_las, shp_fp_rfl,
                 imu_data, known_rfl, json_dir,
                 results_dir, ice_dir, in_dir):
    '''
    
    Runs normalized_reflectance for each of the flight lines and merges into one raster.

    Reprojects it correctly based on  already solved for DEM.

    Calls AART snow model to estimate SSA (grain size).

    
    '''

    # Run through each flightline and mosaic reflectance
    rfl_fp_grid = normalized_reflectance(snow_tif, cal_las, 
                                         shp_fp_rfl,
                                         imu_data, 
                                         known_rfl, 
                                         json_dir,
                                         results_dir, 
                                         ice_dir, in_dir)
    
    # Prepare inputs needed for SSA raster (vectorized operation)
    # theta_grid = 180 (perfect backscatter relative to sensor)
    rfl_grid = rio.open_rasterio(rfl_fp_grid, masked=True)
    rfl_grid = rfl_grid.rio.reproject_match(snowoff)
    ssa_grid = rfl_grid.copy()

    # Call AART
    ssa_grid = aart_1064(rfl_grid, cosi=1, g=0.85, b=1.6)

    return ssa_grid