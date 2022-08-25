"""
Takes input directory full of .laz files and filters+classifies them to DTM laz and DTM tif.

Usage:
    ice-pipeline.py <in_dir>
"""

import os
from docopt import docopt

def replace_white_spaces(parent, replace = ''):
    print(f'Warning! About to replace all whitespaces with "{replace}"s in {os.path.abspath(parent)}')
    print('Press y to continue...')
    ans = input()
    if ans.lower() == 'y':
        for path, folders, files in os.walk(parent):
            for f in files:
                os.rename(os.path.join(path, f), os.path.join(path, f.replace(' ', replace)))
            for i in range(len(folders)):
                new_name = folders[i].replace(' ', replace)
                os.rename(os.path.join(path, folders[i]), os.path.join(path, new_name))
                folders[i] = new_name
    else:
        print(f'Passing...')

if __name__ == '__main__':
    args = docopt(__doc__)
    in_dir = args.get('<in_dir>')
    # convert to abspath
    in_dir = os.path.abspath(in_dir)
    replace_white_spaces(in_dir)