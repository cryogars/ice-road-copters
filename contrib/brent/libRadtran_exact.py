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
    METHOD #1 : solving alpha using libRadtran

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
