import logging
import os
import subprocess
from glob import glob
from os.path import isdir, join
import shutil

from snow_pc.common import make_dirs

def replace_white_spaces(in_dir, replace = ''):
    """Remove any white space in the point cloud files. 

    Args:
        in_dir (_type_): in_dir directory of the point cloud files.
        replace (str, optional): Character to replace the white space. Defaults to ''.
    """
    assert isdir(in_dir), f'{in_dir} is not a directory'
    response = input(f'Warning! About to replace whitespaces with "{replace}"s in {os.path.abspath(in_dir)} \n Press y to continue...')
    if response.lower() == 'y':
        for path, folders, files in os.walk(in_dir):
            for f in files:
                os.rename(os.path.join(path, f), os.path.join(path, f.replace(' ', replace)))
            for i in range(len(folders)):
                new_name = folders[i].replace(' ', replace)
                os.rename(os.path.join(path, folders[i]), os.path.join(path, new_name))
                folders[i] = new_name
    else:
        print(f'Passing...')


def las2laz(in_dir: str):
    """Convert all LAS files in a directory to LAZ files.

    Args:
        in_dir (str): The directory containing the LAS files to convert.
    """

    assert isdir(in_dir), f'{in_dir} is not a directory'

    # Get a list of all LAS files in the directory
    las_files = [file for file in os.listdir(in_dir) if file.endswith('.las')]

    # Iterate over each LAS file and convert it to LAZ
    for las_file in las_files:
        input_path = os.path.join(in_dir, las_file)
        output_path = os.path.join(in_dir, os.path.splitext(las_file)[0] + '.laz')
        subprocess.run(['pdal', 'translate', input_path, output_path])
        print(f"Converted {input_path} to {output_path}")

        
def merge_laz_files(in_dir, out_fp = 'unaligned_merged.laz'):
    """Merge all LAZ files in a directory into a single LAZ file.
    Args:
        in_dir (_type_): Directory containing the LAZ files to merge.
        out_fp (str, optional): Filename of the merged LAZ file. Defaults to 'unaligned_merged.laz'.
    """
    assert isdir(in_dir), f'{in_dir} is not a directory'
    # out fp to save to
    mosaic_fp = join(in_dir, out_fp)
    # Get a list of all LAZ files in the directory
    laz_files = [file for file in os.listdir(in_dir) if file.endswith('.laz')]
    
    # Build the command to merge all LAZ files into a single file
    command = ['pdal', 'merge']
    for laz_file in laz_files:
        command.append(os.path.join(in_dir, laz_file))
    command.append(mosaic_fp)
    
    # Execute the merge command
    print(f'Running command: {command}')
    subprocess.run(command)
    print(f"Merged {len(laz_files)} LAZ files into {mosaic_fp}")

def prepare_pc(in_dir: str, replace: str = ''):
    """Prepare point cloud data for processing.

    Args:
        in_dir (str): Path to the directory containing the point cloud files.
        replace (str, optional): Character to replace the white space. Defaults to ''.

    Returns:
        str: Path to the merged LAZ file.
    """

    # checks on directory and user update
    assert isdir(in_dir), f'Provided: {in_dir} is not a directory. Provide directory with .laz files.'

    #checks if there is at least one file in the directory
    assert len(glob(join(in_dir, '*'))) > 0, f'No files found in {in_dir}'

    #change to the directory
    print(f"Working in directory: {in_dir}")
    os.chdir(in_dir)

    # set up sub directories
    results_dir = make_dirs(in_dir)

    #check and replace white spaces in file paths
    for file in glob(join(in_dir, '*')):
        if ' ' in file:
            print('White spaces found in file paths. Removing...')
            replace_white_spaces(in_dir, replace)
            break

    #check and convert all LAS files to LAZ
    for file in glob(join(in_dir, '*')):
        if file.endswith('.las'):
            print('LAS files found. Converting to LAZ...')
            las2laz(in_dir)
            break
    
    # mosaic
    # if there is more than 1 laz file, merge them
    if len(glob(join(in_dir, '*.laz'))) > 1:
        print('Merging LAZ files...')
        mosaic_fp = os.path.join(results_dir, 'unfiltered_merge.laz')
        merge_laz_files(in_dir, out_fp= mosaic_fp)
        if os.path.exists(mosaic_fp):
            return mosaic_fp
        else:
            print(f"Error: Mosaic file not created")
    #if there is 1 laz file, copy it to the results directory named unfiltered.laz and return the path
    else:
        laz_fp = glob(join(in_dir, '*.laz'))[0]
        shutil.copy(laz_fp, join(results_dir, 'unfiltered.laz'))
        return join(results_dir, 'unfiltered.laz')