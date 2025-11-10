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
import logging
from datetime import datetime
from pathlib import Path


# ----------------------------------------------------------------------
# Logging setup
# ----------------------------------------------------------------------
log = logging.getLogger("ini_template_generator")
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

TEMPLATE = """# ==============================================================
# Ice Road .INI file Configuration Template
# Generated on {date}
# ==============================================================================
# Notes:
#  - All paths can be absolute or relative to this config file.
#  - Leave optional fields blank if not applicable.
#  - Fields marked <required: ...> must be provided before running the pipeline.
# ==============================================================================

[general]
# Required: directory containing .laz or .las LiDAR files
input_dir = <required: path to input LiDAR directory>
# Optional: enable verbose logging (true/false)
debug = false

[dem]
# Optional: user-supplied DEM path. Leave blank to use internal DEM.
user_dem =
# Required if user_dem is provided: true if geoid-based, false if ellipsoid
is_geoid = 

[alignment]
# Required: Shapefile (.shp) to align point clouds with (e.g. ice road centerline)
shapefile = <required: path to alignment shapefile>
# Optional: buffer width in meters (total width, not radius)
buffer_meters = 3.0
# Required: path to Ames Stereo Pipeline (ASP) binaries
asp_dir = <required: path to ASP /bin directory>

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


def write_template_config(filename: Path, force: bool = False) -> None:

    """
    Write the INI template to the target filename.

    - Prompts before overwrite unless `force=True`.
    """

    log.info(f"Preparing to create INI template at: {filename}")
    filename.parent.mkdir(parents=True, exist_ok=True)

    if filename.exists() and not force:
        log.warning(f"Config file already exists at:\n   {filename}")
        ans = input("Overwrite? [y/N]: ").strip().lower()
        if ans != "y":
            log.info("Aborted — file not overwritten.")
            sys.exit(0)

    filename.write_text(
        TEMPLATE.format(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        encoding="utf-8"
    )
    log.info(f"✓ INI configuration template created at:\n   {filename}")
    print("\nNext steps:")
    print(f"  Edit '{filename.name}' and replace all '<required: ...>' fields and placeholders")

def main(argv: list[str]) -> int:
    if len(argv) < 2 or "-h" in argv or "--help" in argv:
        print(__doc__)
        return 0 if len(argv) >= 2 else 1

    # Parse args (very lightweight)
    force = "--force" in argv
    # First non-flag argument is the output path. Safe because we expect a single non-flag arg
    non_flag_args = [a for a in argv[1:] if not a.startswith("-")]

    if len(non_flag_args) == 0:
        log.error("Error: Missing output filename.\n")
        print(__doc__)
        return 1
    
    if len(non_flag_args) > 1:
        log.error(f"Error: Too many arguments: {non_flag_args}")
        log.error("Expected exactly one output filename.\n")
        print(__doc__)
        return 1

    out_path = Path(non_flag_args[0]).expanduser().resolve()
    write_template_config(out_path, force=force)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))