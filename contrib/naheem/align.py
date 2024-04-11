import os
import json
import subprocess
import pandas as pd
import geopandas as gpd
from os.path import dirname, join, exists, basename, abspath
import rioxarray as rxr 
from rasterstats import point_query, zonal_stats

import shutil
from snow_pc.common import download_dem

def clip_align(laz_fp, buff_shp, align_path, asp_dir):
    """Clip the point cloud to a shapefile.

    Args:
        laz_fp (_type_): _description_
        buff_shp (_type_): _description_
        dem_is_geoid (_type_): _description_
        is_canopy (bool, optional): _description_. Defaults to False.

    Raises:
        Exception: _description_
    """
    #set the working directory
    in_dir = dirname(laz_fp)

    #find the directory that the laz_fp is in


    clipped_pc = join(in_dir, 'clipped_pc.laz')
    json_fp = join(in_dir, 'jsons', 'clip_align.json')


    # Create .json file for PDAL clip
    json_pipeline = {
        "pipeline": [
            laz_fp,
            {
                "type":"filters.overlay",
                "dimension":"Classification",
                "datasource":buff_shp,
                "layer":"buffered_area",
                "column":"CLS"
            },
            {
                "type":"filters.range",
                "limits":"Classification[22:22]"
            },
            clipped_pc
        ]
    }
    with open(json_fp,'w') as outfile:
        json.dump(json_pipeline, outfile, indent = 2)

    subprocess.run(['pdal', 'pipeline', json_fp])               

    # Check to see if output clipped point cloud was created
    if not exists(clipped_pc):
        raise Exception('Output point cloud not created')
    
    # set the dem file path
    dem_fp = join(in_dir, 'dem.tif')
    ref_dem = dem_fp

    #There is need to convert the dem to ellipsoid if it is in geoid


    #call asp pc_align function on road and DEM and output translation/rotation matrix
    align_pc = join(in_dir,'pc-align', basename(align_path)) #set the align files name format
    pc_align_func = join(asp_dir, 'pc_align') #set the path to the pc_align function
    subprocess.run([pc_align_func, '--max-displacement', '5', '--highest-accuracy', ref_dem, clipped_pc, '-o', align_pc]) #run the pc_align function

    # Apply transformation matrix to the entire laz and output points
    initial_tansform = align_pc +  '-transform.txt' #set the transform files name format
    transform_pc = join(in_dir,'pc-transform', basename(align_path))
    subprocess.run([pc_align_func, '--max-displacement', '-1', '--num-iterations', '0', '--initial-transform', 
                    initial_tansform, '--save-transformed-source-points', ref_dem, laz_fp,'-o', transform_pc])
    #print the command that was run
    print([pc_align_func, '--max-displacement', '-1', '--num-iterations', '0', '--initial-transform', 
                    initial_tansform, '--save-transformed-source-points', ref_dem, laz_fp,'-o', transform_pc])

    # Grid the output to a 0.5 meter tif (NOTE: this needs to be changed to 1m if using py3dep)
    transform_laz = transform_pc + '-transform.laz'
    point2dem_func = join(asp_dir, 'point2dem')
    subprocess.run([point2dem_func, transform_laz,'--dem-spacing', '0.5', '--search-radius-factor', '2', '-o', align_path])

    return align_path + '-DEM.tif'

def laz_align(laz_fp, align_file, asp_dir, user_dem = ''):
    """Clip the point cloud to a shapefile.

    Args:
        laz_fp (_type_): _description_
        buff_shp (_type_): _description_
        dem_is_geoid (_type_): _description_
        is_canopy (bool, optional): _description_. Defaults to False.

    Raises:
        Exception: _description_
    """
    #set the working directory
    in_dir = dirname(laz_fp)

    # set the dem file path
    dem_fp = join(in_dir, 'dem.tif')

    #download dem using download_dem() if user_dem is not provided
    if user_dem == '':
        dem_fp, crs, project = download_dem(laz_fp, dem_fp= dem_fp)
    else:
        shutil.copy(user_dem, dem_fp) #if user_dem is provided, copy the user_dem to dem_fp

    #if align file is a shapefile
    if align_file.endswith('.shp'):
        buffer_width = 3
        #create a buffer around the shapefile to clip the point cloud
        gdf = gpd.read_file(align_file)
        gdf['geometry'] = gdf.geometry.buffer(buffer_width / 2) #The buffer_width is the entire width. So, must divide by 2 here to get the right distance from centerline.
        gdf['CLS'] = 22 # Create a new attribute to be used for PDAL clip/overlay
        buff_shp = join(in_dir, 'buffered_area.shp')
        gdf.to_file(buff_shp)

        #remove .tif of the laz_fp path and add -align to the end
        align_path = laz_fp.replace('.laz', '-align')

        #set asp_dir
        if basename(asp_dir) != 'bin':
            asp_dir = join(asp_dir, 'bin')

        align_tif = clip_align(laz_fp, buff_shp=buff_shp, align_path= align_path,  asp_dir = asp_dir)

    #elif the file ends with with csv or excel
    elif align_file.endswith('.csv'):
        dem = rxr.open_rasterio(dem_fp)
        #read the csv file
        cal_data = pd.read_csv(align_file)
        #convert the csv file to a geodataframe
        gdf = gpd.GeoDataFrame(cal_data, geometry=gpd.points_from_xy(cal_data['lon'], cal_data['lat']), crs = "EPSG:4326")
        lidar_crs = dem.rio.crs.to_string()
        gdf = gdf.to_crs(lidar_crs)
        #sample the raster values at the points
        vals = point_query(gdf.geometry, dem.squeeze().values, affine = dem.rio.transform(), nodata = -9999)
        # add the values to a new column in the GeoDataFrame
        gdf['dem'] = vals
        #convert DepthCm to float64
        gdf['DepthCm'] = gdf['DepthCm'].astype('float64')
        #convert Deptcm to meters
        gdf['DepthCm'] = gdf['DepthCm']/100
        #create a column SSH in the geodataframe
        gdf['SSH'] = gdf['DepthCm'] + gdf['dem']
        gdf['x'] = gdf.geometry.x
        gdf['y'] = gdf.geometry.y
        #select the x,y and SSH columns
        gdf = gdf[['x', 'y', 'SSH']]
        #drop na
        gdf = gdf.dropna()
        #save to csv
        gdf.to_csv(join(in_dir, 'cal_data.csv'), index = False)
        #remove .tif of the laz_fp path and add -align to the end
        align_path = laz_fp.replace('.laz', '-align')
        #set asp_dir
        if basename(asp_dir) != 'libexec':
            asp_dir = join(asp_dir, 'libexec')
        subprocess.run([join(asp_dir, 'pc_align'), '--max-displacement', '300', '--highest-accuracy', '--datum', 'WGS_1984', '--save-inv-transformed-reference-points', '--save-transformed-source-points', '--csv-format', '1:easting 2: northing 3: height_above_datum', '--csv-proj4', 'EPSG:32611', '--compute-translation-only', laz_fp, join(in_dir, 'cal_data.csv'), '-o', join(in_dir, 'pc_align', basename(align_path))])
        subprocess.run([join(asp_dir, 'point2dem'), join(in_dir, 'pc_align', basename(align_path)) + '-trans_reference.laz', '--dem-spacing', '0.5', '--search-radius-factor', '2', '-o', align_path])
        align_tif = align_path + '-DEM.tif'
    else:
        raise Exception('File type not supported')

    return align_tif






    