"""
Usage:
    ini_template_generator.py <output.ini> [--force] [--no-log-file]
    ini_template_generator.py -h | --help

Description:
    Creates a new configuration template for the Ice Road LiDAR Pipeline.

Arguments:
    <output.ini>
        Target path for the generated INI file. This may be:
          • An absolute path:      /abs/path/run01.ini
          • A relative path:       configs/run01.ini  (relative to current working directory)
          • A bare filename:       run01.ini          (saved in current working directory)

Options:
    --force
        Overwrite the file if it already exists (skips confirmation).
    --no-log-file
        Disable writing to ini_template_generator.log (console output only).
    -h, --help
        Show this help message and exit.

Examples:
    python ini_template_generator.py my_config.ini
    python ini_template_generator.py configs/my_config.ini --force
    python ini_template_generator.py my_config.ini --no-log-file
    python ini_template_generator.py --help
"""

from __future__ import annotations

import sys
import logging
from datetime import datetime
from pathlib import Path


# ----------------------------------------------------------------------
# Template Definition
# ----------------------------------------------------------------------

TEMPLATE = """# ==============================================================================================================================================
# Ice Road .INI file Configuration Template
# Generated on {date}
# ==============================================================================================================================================
# USAGE:
#     • Fill in all fields marked <required>.
#     • Leave optional fields blank unless needed.
#     • Paths may be absolute or relative to this INI file.
#
# FIELD DEFINITIONS:
#     general.input_dir .............. (required): directory containing raw .laz/.las files
#     general.debug .................. (optional: default = false): true/false — enable verbose logging
#
#     dem.user_dem ................... (optional): path to user-specified DEM
#     dem.is_geoid ................... (required IF user_dem is supplied): 
#                                        true  → DEM is geoid-based
#                                        false → DEM is ellipsoid-based
#
#     alignment.shapefile ............ (required): .shp file to align with
#     alignment.buffer_meters ........ (optional: default = 3.0): total width (meters) for transform area
#     alignment.asp_dir .............. (optional; defaults to `./ASP/bin` when not given): 
#                                        directory with ASP binary files
#
#     reflectance.shp_fp_rfl ......... (optional): shapefile to align for reflectance calibration.
#                                        If given, it is assumed you want grain size output.
#                                        Additionally, if this mode is selected, the supplied files must be .LAS with extra bytes included with
#                                        "Intensity as Reflectance" returned by RIEGL.
#     reflectance.imu_data ........... (optional): path to helicopter IMU .CSV or.TXT data used to match data with point cloud using GPS time.
#                                        Column names must include ['Time[s]', 'Easting[m]', 'Northing[m]', 'Height[m]'] 
#     reflectance.cal_las ............ (optional): path to .LAS used for calibration of the apparent reflectance for 1064nm of lidar sensor.
#                                        To avoid confusion, please supply this file in a different directory from <in_dir>.
#     reflectance.known_rfl .......... (optional, float/real): known intrinsic reflectance at 1064 nm for target identified in shp_fp_rfl
#     reflectance.h2o ................ (optional, float): water column vapor in atmosphere (mm)
#     reflectance.aod ................ (optional, float): aerosol optical depth at 550 nm
#
#     smrf.scalar .................... (optional, float): SMRF scalar parameter override
#     smrf.slope ..................... (optional, float): SMRF slope override
#     smrf.threshold ................. (optional, float): SMRF elevation threshold override
#     smrf.window .................... (optional, float): SMRF window size override
# ==============================================================================================================================================

[general]
input_dir = <required>
debug = false

[dem]
user_dem =
is_geoid = <required if user_dem>

[alignment]
shapefile = <required>
buffer_meters = 3.0
asp_dir =

[reflectance]
shp_fp_rfl =
imu_data =
cal_las =
known_rfl =
h2o =
aod =

[smrf]
scalar =
slope =
threshold =
window = 
"""

# ----------------------------------------------------------------------
# Logging setup
# ----------------------------------------------------------------------

def setup_logging(enable_file: bool = True) -> logging.Logger:
    """
    Configure logging for ini_template_generator.

    - Logs to both console and (optionally) a timestamped file in ./logs/
    - File logging can be disabled with enable_file=False
    """
    
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    # Prevent duplicate handlers on repeated runs
    if log.hasHandlers():
        return log

    # Console handler (always active)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    log.addHandler(console_handler)

    if enable_file:
        log_dir = Path.cwd() / "logs"
        log_dir.mkdir(exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"ini_template_generator_{timestamp}.log"

        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        log.addHandler(file_handler)

        log.info(f"Logging to: {log_file}")

    return log

# ----------------------------------------------------------------------
# Core functionality
# ----------------------------------------------------------------------


def write_template_config(
        filename: Path, log: logging.Logger, force: bool = False
) -> None:

    """
    Write the INI template to the target filename.

    - Prompts before overwrite unless `force=True`.
    """

    log.info(f"Preparing to create INI template at:\n   {filename}")
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
    log.info(f"✓ Done!")
    print("\nNext steps:")
    print(f"  Edit '{filename.name}' and replace all '<required: ...>' fields and placeholders")


# ----------------------------------------------------------------------
# CLI Entrypoint
# ----------------------------------------------------------------------

def main(argv: list[str]) -> int:

    """CLI Entrypoint"""

    if len(argv) < 2 or "-h" in argv or "--help" in argv:
        print(__doc__)
        return 0 if len(argv) >= 2 else 1

    # Parse args
    force = "--force" in argv
    no_log = "--no-log-file" in argv
    log = setup_logging(enable_file=not no_log)

    # First non-flag argument = output path
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
    write_template_config(filename=out_path, log=log, force=force)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))