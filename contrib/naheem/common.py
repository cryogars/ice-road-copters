import os
from glob import glob
import pyproj
import laspy
import py3dep
from shapely.geometry import box
from shapely.ops import transform
from rasterio.enums import Resampling
import rioxarray as rxr
from rasterstats import point_query, zonal_stats
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error
from rasterio.crs import CRS
import json
import subprocess


def download_dem(laz_fp, dem_fp, cache_fp ='./cache/aiohttp_cache.sqlite'):
    """Download DEM within the bounds of the las file.

    Args:
        laz_fp (_type_): Path to the las file.
        dem_fp (str, optional): Filename for the downloaded dem. Defaults to 'dem.tif'.
        cache_fp (str, optional): Cache filepath. Defaults to './cache/aiohttp_cache.sqlite'.

    Returns:
        _type_: The filepath to the downloaded DEM, the crs of the las file, and the transform from the las crs to wgs84. 
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)
    
    # read crs of las file
    with laspy.open(laz_fp) as las:
        hdr = las.header
        crs = hdr.parse_crs()
    # log.debug(f"CRS used is {crs}")
    # create transform from wgs84 to las crs
    wgs84 = pyproj.CRS('EPSG:4326')
    project = pyproj.Transformer.from_crs(crs, wgs84 , always_xy=True).transform
    # calculate bounds of las file in wgs84
    utm_bounds = box(hdr.mins[0], hdr.mins[1], hdr.maxs[0], hdr.maxs[1])
    wgs84_bounds = transform(project, utm_bounds)
    # download dem inside bounds
    os.environ["HYRIVER_CACHE_NAME"] = cache_fp
    
    dem_wgs = py3dep.get_map('DEM', wgs84_bounds, resolution=1, crs='EPSG:4326')
    # log.debug(f"DEM bounds: {dem_wgs.rio.bounds()}. Size: {dem_wgs.size}")
    # reproject to las crs and save
    dem_utm = dem_wgs.rio.reproject(crs, resampling = Resampling.cubic_spline)
    dem_utm.rio.to_raster(dem_fp)
    # log.debug(f"Saved to {dem_fp}")
    return dem_fp, crs, project

def make_dirs(in_dir):
    """Create directories for the laz file and the results.

    Args:
        laz_fp (_type_): _description_

    Returns:
        _type_: _description_
    """
    # set up sub directories
    snowpc_dir = os.path.join(in_dir, 'snow-pc')
    os.makedirs(snowpc_dir, exist_ok= True)
    results_dir = os.path.join(snowpc_dir, 'results')
    os.makedirs(results_dir, exist_ok= True)
    return results_dir

def snowdepth_val(lid_path, csv_path, snowdepth_col, lat_col, lon_col, road_shp = '', csv_EPSG=4326, zone_utmcrs=32611, lid_unit="m", probe_unit="cm", use_buffer = 'no'):
    """_summary_

    Args:
        lid_path (_type_): _description_
        csv_path (_type_): _description_
        snowdepth_col (_type_): _description_
        lat_col (_type_): _description_
        lon_col (_type_): _description_
        csv_EPSG (int, optional): _description_. Defaults to 4326.
        zone_utmcrs (int, optional): _description_. Defaults to 32611.
        lid_unit (str, optional): _description_. Defaults to "m".
        probe_unit (str, optional): _description_. Defaults to "cm".
        use_buffer (str, optional): _description_. Defaults to 'no'.

    Returns:
        _type_: _description_
    """
    # read the lidar raster data
    lidar = rxr.open_rasterio(lid_path, masked = True)
    #reproject to crs of the zone
    if lidar.rio.crs.to_string() != "EPSG:" + str(zone_utmcrs):
        lidar = lidar.rio.reproject(CRS.from_string("EPSG:" + str(zone_utmcrs)))
    # read the csv
    df = pd.read_csv(csv_path, usecols=[snowdepth_col, lat_col, lon_col])
    # convert to geodataframe
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df[lon_col], df[lat_col]),
        crs="EPSG:" + str(csv_EPSG),
    )
    # convert the gdf to crs of the zone
    gdf_utm = gdf.to_crs("EPSG:" + str(zone_utmcrs))

    #Extract the LiDAR pixels at pixel or buffered region 
    if use_buffer == "no":
        # sample the snow depth raster values at point locations
        vals = point_query(gdf_utm.geometry, lidar.squeeze().values, affine = lidar.rio.transform(), nodata = -9999)
        # add the values to a new column in the GeoDataFrame
        gdf_utm['lidar'] = vals
    else:
        # buffer the points
        gdf_utm_buffered = gdf_utm.buffer(2)
        # calculate zonal statistics of the mean within the buffer
        vals = zonal_stats(
            gdf_utm_buffered,
            lidar.squeeze().values,
            affine=lidar.rio.transform(),
            nodata=-9999,
            stats="mean",
        )
        # add the values to new columns in the GeoDataFrame
        gdf_utm['lidar'] = [stat["mean"] for stat in vals]

    #convert the unit to m
    if probe_unit == "cm":
        gdf_utm[snowdepth_col] = gdf_utm[snowdepth_col] / 100
    if lid_unit == "cm":
        gdf_utm["lidar"] = gdf_utm["lidar"]/100
    #rename the columns
    gdf_utm.rename(columns={
                       snowdepth_col: 'Probed Snow Depth (m)', 'lidar': 'LiDAR Snow Depth (m)'}, inplace=True)
    #add error column
    gdf_utm['error (cm)'] = (gdf_utm['Probed Snow Depth (m)'] - gdf_utm['LiDAR Snow Depth (m)']) * 100

    if road_shp == '':
        road = '/SNOWDATA/IDALS/misc_data_scripts/3mroadBufferClip/3m_road.shp'
    else:
        road = road_shp
    #read the road shapefile
    road = gpd.read_file(road)
    # convert the gdf to crs of the zone
    road = road.to_crs("EPSG:" + str(zone_utmcrs))
    #create a 5m buffer
    road_10 = road.buffer(2)
    #sample lidar values overlapping with road_10
    values = zonal_stats(
            road_10,
            lidar.squeeze().values,
            affine=lidar.rio.transform(),
            raster_out=True, stats="mean"
        )
    masked_data = values[0]['mini_raster_array']
    lidar_road = masked_data[~masked_data.mask] 

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
    # Drop missing values
    gdf_utm = gdf_utm.dropna(subset=['Probed Snow Depth (m)', 'LiDAR Snow Depth (m)'])

    # Scatter plot of probed depth vs lidar depth
    sns.regplot(x='Probed Snow Depth (m)', y='LiDAR Snow Depth (m)', data=gdf_utm, ax=axs[0], fit_reg=False)
    probe_min, probe_max = gdf_utm['Probed Snow Depth (m)'].min(), gdf_utm['Probed Snow Depth (m)'].max()
    lidar_min, lidar_max = gdf_utm['LiDAR Snow Depth (m)'].min(), gdf_utm['LiDAR Snow Depth (m)'].max()
    min_val = min(probe_min, lidar_min)
    max_val = max(probe_max, lidar_max)
    axs[0].plot([min_val, max_val], [min_val, max_val], linestyle='--', color='black')


    # Calculate correlation (r) and RMSE
    corr = np.corrcoef(gdf_utm['Probed Snow Depth (m)'], gdf_utm['LiDAR Snow Depth (m)'])[0, 1]
    correlation_text = f'r: {corr:.2f}'
    rmse = np.sqrt(mean_squared_error(gdf_utm['Probed Snow Depth (m)'], gdf_utm['LiDAR Snow Depth (m)']))
    rmse_text = f'RMSE: {rmse:.2f}'
    mbe = np.mean(gdf_utm['Probed Snow Depth (m)'] - gdf_utm['LiDAR Snow Depth (m)'])
    mbe_text = f'MBE: {mbe:.2f}'
    mae = np.mean(np.abs(gdf_utm['Probed Snow Depth (m)'] - gdf_utm['LiDAR Snow Depth (m)']))
    mae_text = f'MAE: {mae:.2f}'
    nmad = 1.4826 * np.median(np.abs(gdf_utm['Probed Snow Depth (m)'] - gdf_utm['LiDAR Snow Depth (m)']))
    nmad_text = f'NMAD: {nmad:.2f}'

    # Add the correlation coefficient and RMSE as text to the upper left corner
    combined_text = '\n'.join([correlation_text, rmse_text, mbe_text, mae_text, nmad_text])
    axs[0].text(0.02, 0.98, combined_text, transform=axs[0].transAxes, color='black', bbox=dict(facecolor='white', alpha=0.8), ha='left', va='top')

    # Distribution plot of error
    sns.histplot(gdf_utm['error (cm)'], kde=True, ax=axs[1])
    axs[1].set_xlabel('Error (cm)')
    axs[1].set_ylabel('Frequency')
    plt.tight_layout()


    return gdf_utm, lidar_road

def clip_lidar_with_shapefile(shapefile_path, lidar_input_path, lidar_output_path):
    # Load the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Convert the first geometry in the shapefile to WKT format
    polygon = gdf['geometry'][0].wkt

    # Create a PDAL pipeline to clip the LiDAR file
    pipeline = {
        "pipeline": [
            lidar_input_path,
            {
                "type": "filters.crop",
                "polygon": polygon
            },
            lidar_output_path
        ]
    }

    # Write the pipeline to a JSON file
    with open('pipeline.json', 'w') as f:
        json.dump(pipeline, f)

    # Run the PDAL pipeline
    subprocess.run(['pdal', 'pipeline', 'pipeline.json'])