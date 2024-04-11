"""Main module."""

import os
from os.path import join, basename
import rioxarray as rio

#local imports
from snow_pc.prepare import prepare_pc
from snow_pc.modeling import terrain_models, surface_models
from snow_pc.align import laz_align



def pc2uncorrectedDEM(in_dir, outlas = '', outtif = '', user_dem = '', dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes'):
    """Converts laz files to uncorrected DEM.

    Args:
        in_dir (str): Path to the directory containing the point cloud files.
        user_dem (str, optional): Path to the DEM file. Defaults to ''.

    Returns:
    outtif (str): filepath to output DTM tiff
    outlas (str): filepath to output DTM laz file
    """

    # prepare point cloud
    unfiltered_laz = prepare_pc(in_dir)

    #create uncorrected DTM
    dtm_laz, dtm_tif = terrain_models(unfiltered_laz, outtif= outtif, outlas= outlas, user_dem = user_dem, dem_low = dem_low, dem_high = dem_high, mean_k = mean_k, multiplier = multiplier, lidar_pc = lidar_pc)

    #create uncorrected DSM
    dsm_laz, dsm_tif = surface_models(unfiltered_laz, outtif= outtif, outlas= outlas, user_dem = user_dem, dem_low = dem_low, dem_high = dem_high, mean_k = mean_k, multiplier = multiplier, lidar_pc = lidar_pc)

    return dtm_laz, dtm_tif, dsm_laz, dsm_tif


def pc2correctedDEM(in_dir, align_file, asp_dir, user_dem = '', dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes'):
    """Converts laz files to corrected DEM.

    Args:
        in_dir (str): Path to the directory containing the point cloud files.
        align_shp (str): Path to the shapefile to align the point cloud to.
        user_dem (str, optional): Path to the DEM file. Defaults to ''.

    Returns:
    outtif (str): filepath to output DTM tiff
    outlas (str): filepath to output DTM laz file
    """


    # create uncorrected DEM
    dtm_laz, dtm_tif, dsm_laz, dsm_tif = pc2uncorrectedDEM(in_dir, user_dem = user_dem, dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes')


    # align the point cloud
    dtm_align_tif = laz_align(dtm_laz, align_file = align_file, asp_dir= asp_dir, user_dem = user_dem)
    dsm_align_tif = laz_align(dsm_laz, align_file = align_file, asp_dir= asp_dir, user_dem = user_dem)

    return dtm_align_tif, dsm_align_tif

def pc2snow(in_dir, align_file, asp_dir, user_dem = '', dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes'):
    """Converts laz files to snow depth and canopy height.

    Args:
        in_dir (str): Path to the directory containing the point cloud files.
        align_shp (str): Path to the shapefile to align the point cloud to.
        user_dem (str, optional): Path to the DEM file. Defaults to ''.

    Returns:
    outtif (str): filepath to output DTM tiff
    outlas (str): filepath to output DTM laz file
    """

    # create corrected DEM
    dtm_align_tif, dsm_align_tif = pc2correctedDEM(in_dir, align_file, asp_dir, user_dem = user_dem, dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes')
    
    #set dem_fp
    in_dir = os.path.dirname(dtm_align_tif)
    ref_dem_path = join(in_dir, 'dem.tif')

    #create snow depth
    snow_depth_path = join(in_dir, f'{basename(in_dir)}-snowdepth.tif')
    snowoff = rio.open_rasterio(ref_dem_path, masked=True)
    snowon = rio.open_rasterio(dtm_align_tif, masked=True) 
    snowon_matched = snowon.rio.reproject_match(snowoff)
    snowdepth = snowon_matched - snowoff
    snowdepth.rio.to_raster(snow_depth_path)

    #create canopy height
    canopy_height_path = join(in_dir, f'{basename(in_dir)}-canopyheight.tif')
    canopy_height = rio.open_rasterio(dsm_align_tif, masked=True) - snowoff
    canopy_height.rio.to_raster(canopy_height_path)

    return snow_depth_path, canopy_height_path

