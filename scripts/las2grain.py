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

# NOTE TO SELF:
# Build a LUT of libRadtran runs to facilitate use of alpha param, without BSU folks having to run the radiative transfer
# This should be pretty straight forward. And the altitude distance should be roughly the same for flights.
# I am sure i can get an r2 pretty high using an interpolation method like this..

# ALSO, should consider switching to AOD-550 as more commmon with remote sensing products.
# This may allow us to use a spatially varying AOD & H20 params if/when available.


def calc_transmittance(altitude_km,
                       file_name,
                       path_to_libradtran_bin='/Users/brent/Documents/Albedo/libRadtran-2.0.4/bin', 
                       lrt_dir='/Users/brent/Code/ice-road-copters/test', 
                       path_to_libradtran_base='/Users/brent/Documents/Albedo/libRadtran-2.0.4/',
                       atmos='mw',
                       h=4.79 ,
                       aod=0.073):
    
    '''

    NOTE: h - H20 vapor in mm from AERONET
          aod - is aerosol optical depth at 870 nm from AERONET 
          altitude_km is the altitude of the AERONET station in Merdian, ID
          atmos - assuming mid latitude winter

          https://aeronet.gsfc.nasa.gov/

    Both are from the Meridian_DEQ site. For now I have been updating h and aod variables for each date.
    

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


def grain_pipeline(cal_las, shp_fp_rfl,imu_data, known_rfl,
                   results_dir, ice_dir,in_dir, snow_tif, 
                   snowdepth_fp, canopy_fp):
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
    alpha = estimate_attenuation()  
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
    This is a polynomial approximation that behaves quite well when compared to AART.
    R2 = 1.0 ;mean absolute error = 0.19 μm. 

    '''

    grain_grid = 18008 * (rfl_grid**6)  - 82840 * (rfl_grid**5)  + 161006 * (rfl_grid**4) - 171899 * (rfl_grid**3)  + 109726 * (rfl_grid**2)- 41764 * (rfl_grid**1)  + 7774.1  

    return grain_grid
