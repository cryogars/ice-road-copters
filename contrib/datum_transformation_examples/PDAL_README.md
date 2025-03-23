# PDAL examples to transform datum

Useful for geoid transformations in point cloud data.

See the following for documentation: https://pdal.io/en/stable/stages/filters.reprojection.html .

---------------------------------------------------------------------------------------------------------------

Note: that you may need to consider variable parameters for Helmert transform if transforming between ellipsoids (e.g., WGS to NAD83).

---------------------------------------------------------------------------------------------------------------

For obtaining a geoid grid raster for transformation. If you have installed Ames Stereo Pipeline they provide some. See /ASP/share/geoids. Current options include egm96, egm2008, and navd88. 

---------------------------------------------------------------------------------------------------------------

PDAL can be ran by calling the JSON config file by:

`pdal pipeline -i path/to/your/config.json -v 8`

---------------------------------------------------------------------------------------------------------------