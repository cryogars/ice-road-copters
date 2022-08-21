import py3dep
import rioxarray as rio


# bounding box of the form (west, south, east, north).
geom = (607140, 4865656, 607904, 4866835)

dem = py3dep.get_map('DEM',geom, resolution=1, geo_crs='EPSG:32611')

dem.rio.to_raster("dem.tif")

