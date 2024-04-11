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


def pc2correctedDEM(in_dir, align_file, asp_dir, user_dem = ''):
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
    dtm_laz, dtm_tif, dsm_laz, dsm_tif = pc2uncorrectedDEM(in_dir, user_dem = user_dem)


    # align the point cloud
    dtm_align_tif = laz_align(dtm_laz, align_file = align_file, asp_dir= asp_dir, user_dem = user_dem)
    dsm_align_tif = laz_align(dsm_laz, align_file = align_file, asp_dir= asp_dir, user_dem = user_dem)

    return dtm_align_tif, dsm_align_tif

def pc2snow(in_dir, align_file, asp_dir, user_dem = ''):
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
    dtm_align_tif, dsm_align_tif = pc2correctedDEM(in_dir, align_file, asp_dir, user_dem = user_dem)
    
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


# class Map(ipyleaflet.Map):
#     """Custom map class that inherits from ipyleaflet.Map.
#     """
#     def __init__(self, *args, **kwargs):

#         if "scroll_wheel_zoom" not in kwargs:
#             kwargs["scroll_wheel_zoom"] = True
#         super().__init__(*args, **kwargs)

#         if "layers_control" not in kwargs:
#             kwargs["layers_control"] = True

#         if kwargs["layers_control"]:
#             self.add_LayerControl()

#         if "fullscreen_control" not in kwargs:
#             kwargs["fullscreen_control"] = True

#         if kwargs["fullscreen_control"]:
#             self.add_fullscreen_control()            

#     def add_search_control(self, position = "topleft", **kwargs):
#         """Add a search control to the map.

#         Args:
#             position (str, optional): Position of the search control. Defaults to "topleft".

#         Returns:
#             _type_: SearchControl object.
#         """
#         if "url" not in kwargs:
#             kwargs["url"] = "https://nominatim.openstreetmap.org/search?format=json&q={s}"
#         search = ipyleaflet.SearchControl(position = position, **kwargs)
#         self.add_control(search)
#         return search

#     def add_LayerControl(self, position = "topright"):
#         """Add a layer control to the map.

#         Args:
#             position (str, optional): Position of the layer control. Defaults to "topright".

#         Returns:
#             _type_: LayerControl object.
#         """
#         layer_control = ipyleaflet.LayersControl(position = position)
#         self.add_control(layer_control)
#         return layer_control
    
#     def add_fullscreen_control(self, position = "topright"):
#         """Add a fullscreen control to the map.

#         Args:
#             position (str, optional): Position of the fullscreen control. Defaults to "topright".

#         Returns:
#             _type_: FullscreenControl object.
#         """
#         fullscreen = ipyleaflet.FullScreenControl(position = position)
#         self.add_control(fullscreen)
#         return fullscreen
    
#     def add_tile_layer(self, url, name, **kwargs):
#         """Add a tile layer to the map.

#         Args:
#             url (str): URL of the tile layer.

#         Returns:
#             _type_: TileLayer object.
#         """
#         tile_layer = ipyleaflet.TileLayer(url = url, name=name, **kwargs)
#         self.add_layer(tile_layer)
#         return tile_layer
    
#     def add_basemap(self, basemap, **kwargs):
#         """Add a basemap to the map.

#         Args:
#             basemap (_type_): A string representing the basemap to add.

#         Raises:
#             ValueError: If the basemap is not recognized.
#         """
#         import xyzservices.providers as xyz
#         if basemap.lower() == "openstreetmap":
#             url = "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
#             self.add_tile_layer(url, name = basemap,**kwargs)
#         elif basemap.lower() == "stamen terrain":
#             url = "https://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png"
#             self.add_tile_layer(url, name = basemap,**kwargs)
#         elif basemap.lower() == "opentopomap":
#             url = "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png"
#             self.add_tile_layer(url, name = basemap,**kwargs)
#         elif basemap.lower() == "satellite":
#             url = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
#             self.add_tile_layer(url, name = basemap,**kwargs)

#         else:
#             try:
#                 basemap = eval(f"xyz.{basemap}")
#                 url = basemap.build_url()
#                 name = basemap["name"]
#                 attribute = basemap["attribution"]
#                 print(url, name)
#                 self.add_tile_layer(url, name, attribution = attribute, **kwargs)
#             except:
#                 raise ValueError(f"Basemap {basemap} not recognized.")

#     def add_geojson(self, data, name = "geojson", **kwargs):
#         """Add a GeoJSON layer to the map.

#         Args:
#             data (_type_): A GeoJSON object.

#         Returns:
#             _type_: GeoJSON object.
#         """

#         if isinstance(data, str):
#             import json
#             with open(data, 'r') as f:
#                 data = json.load(f)
#         geojson = ipyleaflet.GeoJSON(data = data, name = name, **kwargs)
#         self.add_layer(geojson)
#         return geojson
    
#     def add_shp(self, data, name = "shapefile", **kwargs):
        # """Add a shapefile to the map.

        # Args:
        #     data (_type_): A shapefile object.

        # Returns:
        #     _type_: GeoData object.
        # """
        # gdf = gpd.read_file(data)
        # geojson = gdf.__geo_interface__
        # self.add_geojson(geojson, name = name, **kwargs)
        
        
