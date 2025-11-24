"""
cli_utils.py
============

Utility functions used across the ice-road-copters processing pipeline.
"""

import os
import logging
from os.path import join
from typing import Optional
from datetime import datetime


def setup_logger(
        log_dir: str, prefix: str, debug: bool = False
)-> logging.Logger:

    """
    Create and configure a logger with both file and console handlers.

    Parameters
    ----------
    log_dir: Directory where log files will be written
    prefix: Prefix for generated log filenames
    debug: If True, enables DEBUG logging level

    Returns
    -------
    A Configured logger instance
    """

    os.makedirs(log_dir, exist_ok=True)
    logfile = join(log_dir, f"{prefix}-{datetime.now().strftime('%Y%m%d-%H%M%S')}.log")

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    formatter = logging.Formatter("(ice-road-copters %(name)s %(levelname)s) %(message)s")

    # File handler
    fh = logging.FileHandler(logfile)
    fh.setFormatter(formatter)

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    # Add handlers
    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.propagate = False

    logger.info(f"Logging initialized at: {logfile}")
    return logger

def replace_white_spaces(
        path: str, replace: str = "", log: Optional[logging.Logger] = None
) -> None:
    
    """
    Recursively rename files and folders under `path`, replacing whitespace
    with a specified replacement character.

    Parameters
    ----------
    path: Path to directory to sanitize
    replace: Replacement character for whitespace (default "_")
    log: Logger to record actions; falls back to print() if None

    Returns
    -------
    Function performs renaming in place and prints/logs actions.
    """

    def _log(msg: str, level: str = "info") -> None:
        if log:
            getattr(log, level)(msg)
        else:
            print(msg)

    _log(
        f"Replacing all whitespaces with '{replace}' inside {os.path.abspath(path)}",
        level="warning"
    )

    for root, folders, files in os.walk(path):
        # rename files
        for f in files:
            new_name = f.replace(" ", replace)
            if new_name != f:
                os.rename(os.path.join(root, f), os.path.join(root, new_name))
                _log(f"Renamed file: {f} -> {new_name}", "debug")

        # rename directories
        for i, folder in enumerate(folders):
            new_folder = folder.replace(" ", replace)
            if new_folder != folder:
                os.rename(os.path.join(root, folder), os.path.join(root, new_folder))
                folders[i] = new_folder
                _log(f"Renamed folder: {folder} -> {new_folder}", "debug")

    _log("Whitespace replacement complete.")
