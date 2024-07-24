def create_libradtran_LUT(
                       path_to_libradtran_bin='/Users/brent/Documents/Albedo/libRadtran-2.0.4/bin', 
                       lrt_dir='/Users/brent/Code/ice-road-copters/test', 
                       path_to_libradtran_base='/Users/brent/Documents/Albedo/libRadtran-2.0.4/',
                       atmos='mw'):

    '''
    Creating LUT for method #2 above (estimate_attenuation_LUT)

    '''

    h20_range = [1, 3, 7, 10, 25, 50]
    a550_range = [0.01, 0.1, 0.2,  0.6, 1.0]
    print(f'starting {len(h20_range) * len(a550_range)} computations for LUT.')
    # Run libRadtran for set of altitude
    altitude_list = np.arange(0, 3.01, 0.25) # 12 runs

    alphas = []
    aods = []
    h20s = []

    for h in h20_range:
        for aod in a550_range:
            for alt in altitude_list:
                # Run here manually for irrad
                file_name = f'run_{alt}_{aod}_{h}'
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
                    f.write(f'altitude {alt}\n')  # altitude  
                    f.write(f'aerosol_default\n')  # 
                    f.write(f'aerosol_species_file continental_average\n')  # 
                    f.write(f'aerosol_set_tau_at_wvl 870 {aod}\n')  #    
                    f.write(f'output_quantity transmittance\n')  #outputs
                    f.write(f'output_user lambda eglo\n')  #outputs  
                    f.write('quiet')

                cmd = f'{path_to_libradtran_bin}/uvspec < {fname}.INP > {fname}.out'
                subprocess.run(cmd, shell=True, cwd=path_to_libradtran_bin)


            # Interp data at 1064 nm
            t_data = []
            my_wavelengths = np.arange(1060, 1070.1, 1)
            for alt in altitude_list: 
                df_t = pd.read_csv(f'{lrt_dir}/run_{alt}_{aod}_{h}.out', delim_whitespace=True, header=None)
                df_t.columns = ['Wavelength', 'eglo']
                # Compute t_up (upward transmittance)
                fun_t = interpolate.interp1d(df_t['Wavelength'], df_t['eglo'], kind='slinear')
                t_1064 = fun_t(my_wavelengths)[4]
                t_data.append([alt, np.log(t_1064)])
            t_array = np.array(t_data)

            # Estimate slope. This is attenuation [km-1]
            alpha ,_ = np.polyfit(t_array[:,0], t_array[:,1], 1)

            alphas.append(alpha)
            h20s.append(h)
            aods.append(aod)
            print(f'{h},{aod} computed.')

    df = pd.DataFrame(data=h20s, columns=['h'])
    df['aod'] = aods
    df['alpha'] = alphas
    df.to_csv('/Users/brent/Code/ice-road-copters/test/libRadtran_LUT.csv')
    
    return 

#create_libradtran_LUT(path_to_libradtran_bin='/Users/brent/Documents/Albedo/libRadtran-2.0.4/bin', 
#                       lrt_dir='/Users/brent/Code/ice-road-copters/test', 
#                       path_to_libradtran_base='/Users/brent/Documents/Albedo/libRadtran-2.0.4/',
#                       atmos='mw')


