"""
Usage:
    ini_template_generator.py <output.ini> [--force]
    ini_template_generator.py -h | --help

Description:
    Creates a new configuration template for the Ice Road LiDAR Pipeline.

Arguments:
    <output.ini>
        Target path for the generated INI file. This may be:
          • An absolute path:      /abs/path/config.ini
          • A relative path:       configs/run01.ini  (relative to current working directory)
          • A bare filename:       run01.ini          (saved in current working directory)

Options:
    --force
        Overwrite the file if it already exists (skips confirmation).
    -h, --help
        Show this help message and exit.

Examples:
    python ini_template_generator.py my_config.ini
    python ini_template_generator.py configs/ice_road_run01.ini --force
    python ini_template_generator.py --help
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

TEMPLATE = """# ==============================================================
# Ice Road .INI file Configuration Template
# Generated on {date}
# ==============================================================
# Notes:
#  - All paths can be absolute or relative to this config file.
#  - Leave optional fields blank if not applicable.
# ==============================================================

[general]
# Required: directory containing .laz or .las LiDAR files
input_dir = ./input_data
# Optional: enable verbose logging (true/false)
debug = false

[dem]
# Optional: user-supplied DEM path. Leave blank to use internal DEM.
user_dem =
# Required if user_dem is provided: true if geoid-based, false if ellipsoid
is_geoid = 

[alignment]
# Shapefile (.shp) to align point clouds with (e.g. ice road centerline)
shapefile = ./shapefiles/road_centerline.shp
# Buffer width in meters (total width, not radius)
buffer_meters = 3.0
# Path to Ames Stereo Pipeline (ASP) binaries
asp_dir = /usr/local/ASP/bin

[reflectance]
# Optional inputs for reflectance or snow grain size analysis
shp_fp_rfl =
imu_data =
cal_las =
known_rfl =
h2o =
aod =

[smrf]
# Optional SMRF (Simple Morphological Filter) overrides
scalar =
slope =
threshold =
window =
"""

if __name__ == "__main__":
    print(TEMPLATE)