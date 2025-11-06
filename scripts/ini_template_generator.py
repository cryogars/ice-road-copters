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
        Overwrite the file if it already exists (skip confirmation).
    -h, --help
        Show this help message and exit.

Examples:
    python ini_template_generator.py my_config.ini
    python ini_template_generator.py configs/ice_road_run01.ini --force
    python ini_template_generator.py --help
"""
