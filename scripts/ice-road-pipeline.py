"""
Takes input directory full of .laz files and filters+classifies them to DTM laz and DTM tif.

Usage:
    ice-pipeline.py <in_dir> [-d debug] [-a asp_dir]

Options:
    -d debug      turns on debugging logging  [default: True]
    -a asp_dir    Directory with ASP binary files [default: False]
"""

from docopt import docopt
from glob import glob
from os.path import abspath, join, basename
from laz2dem import iceroad_logging, las2uncorrectedDEM
from laz_align import laz_align
from datetime import datetime
import logging
import sys
import os

if __name__ == '__main__':
    start_time = datetime.now()
    # get command line args
    args = docopt(__doc__)
    debug = args.get('-d')
    asp_dir = args.get('-a')
    in_dir = args.get('<in_dir>')
    # convert to abspath
    in_dir = abspath(in_dir)
    # setup logging
    log_dir = join(in_dir, 'logs')
    os.makedirs(log_dir, exist_ok= True)
    log_prefix = 'ice-road-pipeline'
    old_logs = glob(join(log_dir, f'{log_prefix}*.log'))
    if old_logs:
        vnum = max([int(basename(i).split('.')[0].split('-')[-1].replace('r','')) for i in sorted(old_logs)]) + 1
    else:
        vnum = 1

    logging.basicConfig(level=logging.INFO,
    format=f"(ice-road-copters {__name__} %(levelname)s) %(message)s",
    handlers=[
        logging.FileHandler(join(log_dir, f'{log_prefix}-r{vnum}.log')),
        logging.StreamHandler(sys.stdout)]
    )
    log = logging.getLogger(__name__)
    if debug:
        log.setLevel(logging.DEBUG)

    # run main function
    outtif, outlas = las2uncorrectedDEM(in_dir, debug, log)
    if asp_dir == 'False':
        asp_dir = False
    aligned_tif = laz_align(join(in_dir, 'results'), asp_dir = asp_dir)

    end_time = datetime.now()
    log.info(f"Completed! Run Time: {end_time - start_time}")
