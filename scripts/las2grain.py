from glob import glob
from os.path import join, basename
import rioxarray as rio
import rasterio
import numpy as np
import pandas as pd
from scipy import interpolate
import logging
import re
import json
import subprocess
from subprocess import Popen, PIPE
import os

# local imports
from laz2dem import cl_call


log = logging.getLogger(__name__)


def calc_transmittance(altitude_km,
                       file_name,
                       path_to_libradtran_bin='/Users/brent/Documents/Albedo/libRadtran-2.0.4/bin', 
                       lrt_dir='/Users/brent/Code/ice-road-copters/test', 
                       path_to_libradtran_base='/Users/brent/Documents/Albedo/libRadtran-2.0.4/',
                       atmos='mw',
                       h=5.93 ,
                       aod=0.095):
    
    '''

    NOTE: h - H20 vapor in mm from AERONET
          aod - is aerosol optical depth at 870 nm from AERONET 
          altitude_km is the altitude of the AERONET station in Merdian, ID
          atmos - assuming mid latitude winter

          https://aeronet.gsfc.nasa.gov/

    Both are from the Meridian_DEQ site for Feb 9,2023. These values can be changes as needed.
    For now they are hard coded here in the function.

    On the physics, the attenuation is relatively small using our methodology, and only really shows a small 1% 
    change when the range (or beam distance) varies largely from the calibration range at Eagle.

    This function is used in a for loop to solve for extinction coef based on change in altitude and transmittance.




    '''
    

    # Run here manually for irrad
    fname = f'{lrt_dir}/{file_name}'
    with open(f'{fname}.INP', 'w') as f:
        f.write('source solar\n')  # extraterrestrial spectrum
        f.write('wavelength 600 1100\n')  # set range for lambda
        f.write(f'atmosphere_file {path_to_libradtran_base}/data/atmmod/afgl{atmos}.dat\n')
        f.write(f'albedo {0.0}\n')  # 
        f.write(f'umu 1 \n') # Cosine of the view zenith angle
        f.write(f'phi 180 \n') # VAA
        f.write('rte_solver disort\n')  # set disort
        f.write('pseudospherical\n')# computed with spherical instead of plane parallel
        f.write(f'mol_modify O3 300 DU\n')  #  
        f.write(f'mol_abs_param reptran coarse\n')  #  
        f.write(f'mol_modify H2O {h} MM\n')  #  
        f.write(f'crs_model rayleigh bodhaine \n')  # 
        f.write(f'zout sur\n')  # 
        f.write(f'altitude {altitude_km}\n')  # altitude  
        f.write(f'aerosol_default\n')  # 
        f.write(f'aerosol_species_file continental_average\n')  # 
        f.write(f'aerosol_set_tau_at_wvl 870 {aod}\n')  #    
        f.write(f'output_quantity transmittance\n')  #outputs
        f.write(f'output_user lambda eglo\n')  #outputs  
        f.write('quiet')

    cmd = f'{path_to_libradtran_bin}/uvspec < {fname}.INP > {fname}.out'
    subprocess.run(cmd, shell=True, cwd=path_to_libradtran_bin)
    
    return




def estimate_attenuation(lrt_dir='/Users/brent/Code/ice-road-copters/test'):

    '''
    This is a simple approach for solving for attenuation [km-1] where
    the true altitude of Merdian AERONET is 0.808 km , and looping nearby to interpolate change.

    '''

    # Run libRadtran for set of altitude
    altitude_list = np.arange(0, 3.01, 0.25) # 12 runs
    for alt in altitude_list:
        calc_transmittance(alt, f'run_{alt}_km')
    
    # Interp data at 1064 nm
    t_data = []
    my_wavelengths = np.arange(1060, 1070.1, 1)
    for alt in altitude_list: 
        df_t = pd.read_csv(f'{lrt_dir}/run_{alt}_km.out', delim_whitespace=True, header=None)
        df_t.columns = ['Wavelength', 'eglo']
        # Compute t_up (upward transmittance)
        fun_t = interpolate.interp1d(df_t['Wavelength'], df_t['eglo'], kind='slinear')
        t_1064 = fun_t(my_wavelengths)[4]
        t_data.append([alt, np.log(t_1064)])
    t_array = np.array(t_data)

    # Estimate slope. This is attenuation [km-1]
    alpha ,_ = np.polyfit(t_array[:,0], t_array[:,1], 1)

    #print(alpha)
    #import matplotlib.pyplot as plt
    #plt.scatter(t_array[:,0], t_array[:,1])
    #plt.xlabel('Altitude [km]')
    #plt.ylabel(r'$ln(\tau)$')
    #plt.show()

    
    return alpha


def normalized_reflectance(snowon_matched, cal_las, shp_fp_rfl,
                           imu_data, known_rfl, json_dir,
                           results_dir, ice_dir, in_dir, snow_tif):
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
    crs = snowon_matched.crs
    pix_size, _ = snowon_matched.res
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

    # Estimate atmosph attenuation
    alpha = estimate_attenuation()      
    log.info(f'extinction coefficient: {alpha}') 

    # Calcs calibration stats from target
    output_csv = f'{ssa_dir}/all-calibration-rfl.csv'
    cl_call(f'Rscript {scripts_dir}/las_ssa_cal.R {cal_las} {crs} {shp_fp_rfl} {n_e_d_shift} \
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
        
        # Set file path to this specific .las in loop
        rfl_fp = f'{ssa_dir}/{las_name}-rfl.las'

        # Getting a translated "LAS" file
        # "LAS" in quotations bc I am hiding the rfl here in "Z"
        # with the intention to do fast IDW in the next step.
        cl_call(f'Rscript {scripts_dir}/las_ssa_prep.R {f} {crs} {ni_fp} {nj_fp} {nk_fp} \
                {rfl_fp} {n_e_d_shift} {road_cal_factor} {imu_data} {str(alpha)}', log)

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

    crs = snowon_matched.crs

    rfl_grid = rio.open_rasterio(rfl_fp_grid, masked=True)
    rfl_grid = rfl_grid.rio.set_crs(crs, inplace=True)
    rfl_grid.rio.to_raster(rfl_fp_grid)

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





def ssa_pipeline(snowon_matched, cal_las, shp_fp_rfl,
                 imu_data, known_rfl, json_dir,
                 results_dir, ice_dir, in_dir, snow_tif):
    '''
    
    Runs normalized_reflectance for each of the flight lines and merges into one raster.

    Reprojects it correctly based on  already solved for DEM.

    Calls AART snow model to estimate SSA (grain size).

    
    '''

    # Run through each flightline and mosaic reflectance
    rfl_fp_grid = normalized_reflectance(snowon_matched, cal_las, 
                                         shp_fp_rfl,
                                         imu_data, 
                                         known_rfl, 
                                         json_dir,
                                         results_dir, 
                                         ice_dir, in_dir,
                                         snow_tif)
    
    rfl_grid = rio.open_rasterio(rfl_fp_grid, masked=True)
    rfl_grid = rfl_grid.rio.reproject_match(snowon_matched)
    ssa_grid = rfl_grid.copy()

    # Call AART
    ssa_grid = aart_1064(rfl_grid, cosi=1, g=0.85, b=1.6)

    return ssa_grid