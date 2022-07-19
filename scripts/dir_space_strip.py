import os
parent_dir = '../test/'

def replace(parent, replace = ''):
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

replace(parent_dir)