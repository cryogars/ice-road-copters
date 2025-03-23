# PDAL examples to transform datum

Useful for geoid transformations in point cloud data.

See the following for documentation: https://pdal.io/en/stable/stages/filters.reprojection.html .

---------------------------------------------------------------------------------------------------------------

Note: that you may need to consider variable parameters for Helmert transform if transforming between ellipsoids (e.g., WGS to NAD83).

---------------------------------------------------------------------------------------------------------------

PDAL can be ran by calling the JSON config file by:

`pdal pipeline -i path/to/your/config.json -v 8`

---------------------------------------------------------------------------------------------------------------