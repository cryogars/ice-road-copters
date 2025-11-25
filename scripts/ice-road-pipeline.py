"""
Takes input directory full of .laz (or.las) files and filters+classifies them to DTM laz and DTM tif.


Usage:
    ice-road-pipeline.py <in_dir> [--shp=<shpfile>] [--dem=<demfile>]
                       [--debug] [--geoid] [--skip-filter] [--force]
                       [--asp-dir=<dir>] [--buffer=<m>]
                       [--smrf-scalar=<x>] [--smrf-slope=<y>]
                       [--smrf-threshold=<z>] [--smrf-window=<w>]
                       [--shp-rfl=<shp>] [--imu=<csv>] [--cal-las=<las>]
                       [--known-rfl=<val>] [--h2o=<val>] [--aod=<val>]

Arguments:
    <in_dir>                  Directory containing .laz or .las files

Required Arguments:
    --shp=<shpfile>           Shapefile to align with (.shp)
    --dem=<demfile>           Path to user specifed DEM (.tif/.tiff)

Flags:
    --debug                   Enable debug logging (default = False)
    --geoid                   Input DEM uses geoid vertical datum (default = False)
    --skip-filter             Skip DEM filtering after ASP alignment (default = False)
    --force                   Force overwrite existing TIFs (uncorrected and aligned) and strip whitespace from <in_dir>

Optional:
    --asp-dir=<dir>           Directory with ASP binary files (default: ./ASP/bin)
    --buffer=<m>              Buffer width in meters [default: 3.0]

SMRF Overrides (optional):
    --smrf-scalar=<x>         SMRF scalar parameter (float)
    --smrf-slope=<y>          SMRF slope parameter (float)
    --smrf-threshold=<z>      SMRF threshold parameter (float)
    --smrf-window=<w>         SMRF window parameter (float)

Grain-size reflectance (optional):
    --shp-rfl=<shp>           Shapefile to align for reflectance calibration. If given, it is assumed you want grain size output.
                              Additionally, if this mode is selected, the supplied files must be .LAS with extra bytes included with
                              "Intensity as Reflectance" returned by RIEGL.
    --imu=<csv>               Path to helicopter IMU .CSV or.TXT data used to match data with point cloud using GPS time.
                              Column names must include ['Time[s]', 'Easting[m]', 'Northing[m]', 'Height[m]']
    --cal-las=<las>           Path to .LAS used for calibration of the apparent reflectance for 1064nm of lidar sensor.
                              To avoid confusion, please supply this file in a different directory from <in_dir>.
    --known-rfl=<val>         Known intrinsic reflectance at 1064nm (float/real) for target identified in shp_fp_rfl
    --h2o=<val>               Water Column Vapor in atmosphere in mm (float)
    --aod=<val>               Aerosol optical depth at 550nm (float)
"""

from docopt import docopt
from schema import Schema, And, Or, Use, SchemaError
from glob import glob
from os.path import abspath, join, basename, isdir
import rioxarray as rio
from datetime import datetime
import logging
import sys
import os

# local imports
from laz2dem import las2uncorrectedDEM
from laz_align import laz_align
from las2grain import grain_pipeline
from utils import setup_logger, replace_white_spaces

SMRF_OPTIONS = [
    ('--smrf-scalar', 'scalar'),
    ('--smrf-slope', 'slope'),
    ('--smrf-threshold', 'threshold'),
    ('--smrf-window', 'window'),
]

def main():

    start_time = datetime.now()

    # ---- docopt parse and in_dir validation ----
    try:
        args = docopt(__doc__)
    except SystemExit:
        # If user requested --help, allow docopt to print full help text normally
        if "--help" in sys.argv or "-h" in sys.argv:
            sys.exit(0)

        print("\nERROR: Missing input directory <in_dir>.\n")
        print("Example:\n  ice-road-pipeline.py /tmp/data --shp=road.shp --dem=ref.tif\n")
        sys.exit(1)

    # ---- Enforce user DEM(-e) and alignment shapefile(-e) ----
    required = ["--shp", "--dem"]
    friendly = {
        "--shp": "alignment shapefile (--shp=<file>)",
        "--dem": "reference DEM (--dem=<file>)",
    }

    missing = [r for r in required if not args.get(r)]
    if missing:
        print(f"\nERROR: Missing required argument: {friendly[missing[0]]}\n")
        sys.exit(1)

    # ---- schema validation ----
    schema = Schema({
        "<in_dir>": And(os.path.exists, os.path.isdir,
                        error="<in_dir> must be an existing directory"),
        "--shp": And(os.path.exists, lambda p: p.lower().endswith('.shp'),
                     error="--shp must be a .shp file"),
        "--dem": And(os.path.exists, lambda p: p.lower().endswith(('.tif', '.tiff')),
                     error="--dem must be a .tif/.tiff DEM"),

        "--buffer": Or(None, Use(float, error="--buffer must be numeric")),
        "--smrf-scalar": Or(None, Use(float, error="--smrf-scalar must be numeric")),
        "--smrf-slope": Or(None, Use(float, error="--smrf-slope must be numeric")),
        "--smrf-threshold": Or(None, Use(float, error="--smrf-threshold must be numeric")),
        "--smrf-window": Or(None, Use(float, error="--smrf-window must be numeric")),

        "--debug": bool,
        "--geoid": bool,
        "--skip-filter": bool,
        "--force": bool,

        "--asp-dir": Or(None, str),
        "--shp-rfl": Or(None, str),
        "--imu": Or(None, str),
        "--cal-las": Or(None, str),
        "--known-rfl": Or(None, Use(float)),
        "--h2o": Or(None, Use(float)),
        "--aod": Or(None, Use(float)),
    })

    try:
        args = schema.validate(args)
    except SchemaError as e:
        msg = str(e)
        if " in " in msg:              # remove verbose dict printing
            msg = msg.split(" in ")[0]
        
        print("\nArgument Error:")
        print(f"  {msg}\n")
        sys.exit(1)
    
    # ---- extract args ----
    in_dir      = abspath(args["<in_dir>"])
    shp_fp      = abspath(args["--shp"])
    user_dem    = abspath(args["--dem"])

    debug       = args["--debug"]
    geoid       = args["--geoid"]
    skip_filter = args["--skip-filter"]
    force       = args["--force"]
    asp_dir     = args["--asp-dir"]

    if asp_dir:
        asp_dir = abspath(asp_dir)
        if basename(asp_dir) != "bin":
            asp_dir = join(asp_dir, "bin")
    else:
        asp_dir = abspath(join("ASP", "bin"))

    buffer_meters  = args.get("--buffer", 3.0)
    smrf_overrides = {key: args[flag] for flag, key in SMRF_OPTIONS if args[flag] is not None}

    # grain args (todo: split grain vs IR pipeline)
    shp_fp_rfl = args["--shp-rfl"]
    imu_data   = args["--imu"]
    cal_las    = args["--cal-las"]
    known_rfl  = args["--known-rfl"]
    h2o        = args["--h2o"]
    aod        = args["--aod"]
    if shp_fp_rfl:
        shp_fp_rfl = abspath(shp_fp_rfl)
    las_extra_byte_format = bool(shp_fp_rfl)

    # ---- directory setup ----
    ice_dir     = join(in_dir, "ice-road")
    results_dir = join(ice_dir, "results")
    json_dir    = join(ice_dir, "jsons")
    log_dir     = join(ice_dir, "logs")

    for d in (ice_dir, results_dir, json_dir, log_dir):
        os.makedirs(d, exist_ok=True)

    # ---- logging ----
    log = setup_logger(log_dir, "ice-road-pipeline", debug)
    log.info("Arguments validated. Starting IRC pipeline.")

    # ---- whitespace cleanup ----
    if any(" " in p for p in glob(join(in_dir, "*"))):
        log.warning(f"Whitespace(s) found in filenames within {in_dir}.")
        if force:
            replace_white_spaces(path=in_dir, log=log)
        else:
            log.warning("Not modifying (rerun with --force to rename).")

    # ---- DEM creation ----
    log.info("Starting laz2uncorrectedDEM...")
    log.info(f'Using in_dir: {in_dir}, user_dem: {user_dem}')
    outtif, outlas, canopy_laz = las2uncorrectedDEM(
        in_dir=in_dir,
        log=log,
        debug=debug,
        las_extra_byte_format=las_extra_byte_format,
        smrf_overrides=smrf_overrides,
        force_overwrite=force,
    )

    # ---- ASP alignment ----
    log.info('Starting ASP laz align')
    log.info(f'Using in_dir: {in_dir}, shapefile: {shp_fp}, ASP dir: {asp_dir}')
    snow_tif, canopy_tif = laz_align(
        in_dir=in_dir,
        align_shp=shp_fp,
        asp_dir=asp_dir,
        log=log,
        input_laz=outlas,
        canopy_laz=canopy_laz,
        dem_is_geoid=geoid,
        buffer_meters=buffer_meters,
        use_dem_filter=not skip_filter,
        user_dem=user_dem,
        force_overwrite=force,
    )

    # cleanup up after ASP a bit
    for fp in os.listdir(ice_dir):
        if fp.endswith(".txt"):
            os.remove(join(ice_dir, fp))
        if fp.endswith("-DEM.tif"):
            os.rename(join(ice_dir, fp), join(ice_dir, fp.replace("-DEM", "")))
    
    snow_tif   = snow_tif.replace("-DEM", "")
    canopy_tif = canopy_tif.replace("-DEM", "")

    # ---- snow depth ----
    # difference two rasters to find snow depth
    ref_dem_path = join(results_dir, 'dem.tif')
    snow_depth_path = join(ice_dir, f'{basename(in_dir)}-snowdepth.tif')

    snowoff = rio.open_rasterio(ref_dem_path, masked=True)
    snowon = rio.open_rasterio(snow_tif, masked=True)
    snowon_matched = snowon.rio.reproject_match(snowoff)
    snowdepth = snowon_matched - snowoff

    # Write reproject-match SnowDepth and SnowOn to disk
    snowdepth.rio.to_raster(snow_depth_path)
    snowon_matched.rio.to_raster(snow_tif)

    # ---- canopy height ----
    # difference two rasters to find canopy height
    canopy_fp    = join(ice_dir, f'{basename(in_dir)}-canopyheight.tif')
    canopy       = rio.open_rasterio(canopy_tif, masked=True)
    matched      = canopy.rio.reproject_match(snowoff)
    canopyheight = matched - snowoff

    # mask snow depth from vegetation (null=0)
    canopyheight = canopyheight.where((canopyheight > snowdepth + 0.1) | (snowdepth.isnull()), other=0)
    canopyheight.rio.to_raster(canopy_fp)

    # ---- grain pipeline ----
    if shp_fp_rfl:
        log.info("Running grain size pipeline...")
        grain_pipeline(
            cal_las, shp_fp_rfl,
            imu_data, known_rfl,
            results_dir, ice_dir,
            in_dir, snow_tif,
            snow_depth_path,
            canopy_fp, h2o, aod
        )
    
    end_time = datetime.now()
    log.info(f"Completed! Pipeline runtime: {end_time - start_time}")

if __name__ == '__main__':
    main()