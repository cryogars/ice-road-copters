from glob import glob
from os.path import join, basename
import rioxarray as rio
import rasterio
import numpy as np
import pandas as pd
from scipy import interpolate
import logging
import re
import subprocess
from subprocess import Popen, PIPE
import os

# local imports
from laz2dem import cl_call


log = logging.getLogger(__name__)



def estimate_attenuation_LUT(h,aod):

    '''
    Solving alpha using LUT

    I built a Lookup-Table using the following range of values:
            h20_range = [1, 3, 7, 10, 25, 50]
            a550_range = [0.01, 0.1, 0.2,  0.6, 1.0]
    
    These 30 simulations are saved within this function, and then are interpolated.

    User inputs:
        h -  water column vapor in mm
        aod - aerosol optical depth at 870 nm

    Outputs:
        alpha - attenuation [km-1]

    '''

    # Alphas solved for previously using libRadtran using mid-lat winter, O3=300 DU, and R[0-3 km] 
    h_arr = np.array([1,  1,  1,  1,  1,  3,  3,  3,  3,  3,  7,  7,  7,  7,  7, 10,
                      10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50])
    aod_arr = np.array([0.01, 0.1,  0.2,  0.6,  1., 0.01, 0.1,  0.2,  0.6,  1.,   
                        0.01, 0.1,  0.2,  0.6, 1., 0.01, 0.1, 0.2,  0.6,  1.,   
                        0.01, 0.1,  0.2,  0.6,  1.,   0.01, 0.1,  0.2,
                        0.6,  1.  ])
    alpha_arr = np.array([0.00332158, 0.00359646, 0.00390284, 0.00516555, 0.00648356, 0.00337639,
                          0.0040808,  0.00487582, 0.00816571, 0.01158536, 0.00360503, 0.00626687,
                          0.00930227, 0.02206811, 0.03556593, 0.0034354,  0.00443985, 0.00558133,
                          0.01035168, 0.01535739, 0.00338004, 0.00341727, 0.00345433, 0.00358356,
                          0.00369771, 0.00343997, 0.00345951, 0.00347671, 0.00352506, 0.00355757])

    #interpolation
    f = interpolate.interp2d(h_arr, aod_arr, alpha_arr, kind = 'cubic')
    alpha = f(h, aod)[0]

    return alpha



def grain_pipeline(cal_las, shp_fp_rfl,imu_data, known_rfl,
                   results_dir, ice_dir,in_dir, snow_tif, 
                   snowdepth_fp, canopy_fp, h2o, aod):
    '''
    This function utilizes lidR package to leverage the use of key features such as translation,
    sample grid data, and apply math operations. With this being an R package, this function requires
    the installation of R and the following packages:

    raster, lidR, rlas, dplyr, readr, terra,sf, data.table

    In summary,
    1 of the n LAS fight data is input, LAS is shifted in X, Y, and Z direction based on ASP pc_align with translation only,
    Heli IMU data (X,Y,Z,t) is used to build vector between surface and sensor --> used to estimate incidence angles, 
    Reflectance is normalized + cleaned, and Optical Grain Size is estimated based on AART.

    '''

    # Check whether Optical Grain Size directory exists, if not, make one
    grain_dir = f'{results_dir}/grain-calc'
    if not os.path.exists(grain_dir):
        os.makedirs(grain_dir)

    # Make sure R is installed
    proc = Popen(["which", "R"],stdout=PIPE,stderr=PIPE)
    exit_code = proc.wait()
    if exit_code != 0:
        raise Exception("Please install R on this system and packages: raster, lidR, rlas, dplyr, readr, terra,sf, data.table.")

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

    # Gdal-warp to 3m spatial resolution for DEM
    dem_spacing = 3
    new_dem_fp = join(ice_dir, f'{basename(in_dir)}-snowon_DEM_{dem_spacing}m.tif')
    os.system(f'gdalwarp -r bilinear -t_srs {crs} \
                -tr {dem_spacing} {dem_spacing} -overwrite {snow_tif} {new_dem_fp} -q')           

    # New DEM
    ref_raster = rasterio.open(new_dem_fp)
    crs2 = ref_raster.crs
    pix_size, _ = ref_raster.res
    crs = str(crs2).split(":",1)[1]
    log.info(f'Using pixel size: {pix_size} m and using CRS: {crs}.') 

    # Compute slope and aspect from snow-on lidar
    slope_fp = join(ice_dir, f'{basename(in_dir)}-snowon_slope.tif')
    aspect_fp = join(ice_dir, f'{basename(in_dir)}-snowon_aspect.tif')
    os.system(f'gdaldem slope -compute_edges {new_dem_fp} {slope_fp} -q')
    os.system(f'gdaldem aspect -compute_edges -zero_for_flat {new_dem_fp} {aspect_fp} -q')
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

    # Estimate atmosph attenuation
    alpha = estimate_attenuation_LUT(h=h2o, 
                                     aod=aod)  
    log.info(f'extinction coefficient: {alpha}') 

    # Calcs calibration stats from target
    output_csv = f'{grain_dir}/all-calibration-rfl.csv'
    cl_call(f'Rscript {scripts_dir}/las_grain_cal.R {cal_las} {crs} {shp_fp_rfl} {n_e_d_shift} \
            {output_csv} {imu_data} {str(pix_size)} {str(alpha)}', log)

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
        
        # Set file path to this specific .tif in loop
        rfl_fp = f'{grain_dir}/{las_name}-rfl.tif'
        rfl_fp_3m = f'{grain_dir}/{las_name}-rfl_3m.tif'
        grain_fp = f'{grain_dir}/{las_name}-grain_size_3m.tif'

        # Getting a RFL tif
        cl_call(f'Rscript {scripts_dir}/las_grain_prep.R {f} {crs} {ni_fp} {nj_fp} {nk_fp} {snow_tif} \
                {rfl_fp} {n_e_d_shift} {road_cal_factor} {imu_data} {str(alpha)} \
                {snowdepth_fp} {canopy_fp}', log)
        
        # Resample RFL tif to 3 m (from 0.5 m)
        os.system(f'gdalwarp -r bilinear -t_srs {crs2} \
            -tr {dem_spacing} {dem_spacing} -overwrite {rfl_fp} {rfl_fp_3m} -q')           

        # Compute Optical Grain Size
        rfl_grid = rio.open_rasterio(rfl_fp_3m, masked=True)
        grain_grid = aart_1064(rfl_grid)
        grain_grid.rio.to_raster(grain_fp)
        

    return


def aart_1064(rfl_grid):


    '''
    Solves Optical Grain Size for AART at 1064 nm wavelength. 
    Assumes pixel is 100% snow and that lidar is pure backscatter.
    g=0.75 and b=1.6 are assumed
    This is an exact solution using algebra of Kokhanovsky et al. (2021).  

    '''
    # Constants
    LAM = 1064*1e-9
    THETA = 180
    COS_SZA = 1
    COS_VZA = 1
    K_EFF = 1.89839284188334E-06
    G = 0.75
    B = 1.6

    # AART
    p = 11.1 * np.exp(-0.087 * THETA) + 1.1 * np.exp(-0.014 * THETA)
    r0 = (1.247 + 1.186 * (COS_SZA + COS_VZA) + 5.157 * COS_SZA * COS_VZA + p) / (4.0 *(COS_SZA + COS_VZA))
    eps =  (9*(1-G)) / (16*B)
    gamma_i = (4 * np.pi * K_EFF) / LAM
    u1 = 0.6*COS_SZA + 1. / 3. + np.sqrt(COS_SZA) / 3.
    u2 = 0.6*COS_VZA + 1. / 3. + np.sqrt(COS_VZA) / 3.
    f = (u1 * u2) / r0

    # Eqn (12) in Wilder et al. (2025)
    grain_grid = (((-np.ln((rfl_grid/r0)**(1/f)))^2) / (2 * eps * gamma_i)) * 1e6

    return grain_grid