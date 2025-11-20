# ASP dem-geoid examples

Useful for geoid transformations in raster data.

See the following for documentation: https://stereopipeline.readthedocs.io/en/latest/tools/dem_geoid.html .

---------------------------------------------------------------------------------------------------------------
### Converting from WGS84 ellipsoid to the EGM2008 geoid
`dem_geoid input-DEM.tif --geoid egm2008`

---------------------------------------------------------------------------------------------------------------

### Converting from EGM2008 geoid to the WGS84 ellipsoid
`dem_geoid input-DEM.tif --geoid egm2008 --reverse-adjustment`

---------------------------------------------------------------------------------------------------------------

### Converting from NAD83 ellipsoid to the NAVD88 geoid
`dem_geoid input-DEM.tif --geoid NAVD88`

---------------------------------------------------------------------------------------------------------------

### Converting from NAVD88 geoid to the NAD83 ellipsoid
`dem_geoid input-DEM.tif --geoid NAVD88 --reverse-adjustment`

---------------------------------------------------------------------------------------------------------------