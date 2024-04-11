import os
from os.path import dirname, join
import json
import subprocess
import shutil
from snow_pc.common import download_dem, make_dirs


#combine the filters into a single function
def terrain_models(laz_fp, outlas = '', outtif = '', user_dem = '', dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes'):
    """Use filters.dem, filters.mongo, filters.elm, filters.outlier, filters.smrf, and filters.range to filter the point cloud for terrain models.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.
        outlas (str, optional): Filepath to save the output las file. Defaults to ''.
        outtif (str, optional): Filepath to save the output tif file. Defaults to ''.
        user_dem (str, optional): Filepath to the dem file. Defaults to ''.
        dem_low (int, optional): _description_. Defaults to 20.
        dem_high (int, optional): _description_. Defaults to 50.
        mean_k (int, optional): _description_. Defaults to 20.
        multiplier (int, optional): _description_. Defaults to 3.

    Returns:
        _type_: Filepath to the terrain model.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    # os.chdir(in_dir)

    #create a filepath for the output las and tif file
    if outlas == '':
        outlas = join(in_dir, 'dtm.laz')
    if outtif == '':
        outtif = join(in_dir, 'dtm.tif')

    #set dem_fp
    dem_fp = join(in_dir, 'dem.tif')
    
    #download dem using download_dem() if user_dem is not provided
    if user_dem == '':
        dem_fp, crs, project = download_dem(laz_fp, dem_fp= dem_fp)
    else:
        shutil.copy(user_dem, dem_fp) #if user_dem is provided, copy the user_dem to dem_fp

    if lidar_pc.lower() == 'yes':
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
                },
                {
                    "type": "filters.dem",
                    "raster": dem_fp,
                    "limits": f"Z[{dem_low}:{dem_high}]"
                },
                {
                    "type": "filters.mongo",\
                    "expression": {"$and": [\
                    {"ReturnNumber": {"$gt": 0}},\
                    {"NumberOfReturns": {"$gt": 0}} ] }
                },
                {
                    "type": "filters.elm"
                },
                {
                    "type": "filters.outlier",\
                    "method": "statistical",\
                    "mean_k": mean_k,\
                    "multiplier": multiplier
                },
                {
                    "type": "filters.smrf",\
                    "ignore": "Classification[7:7], NumberOfReturns[0:0], ReturnNumber[0:0]"
                },
                {
                    "type": "filters.range",
                    "limits": "Classification[2:2]"
                },
                {
                    "type": "writers.las",
                    "filename": outlas,
                    "major_version": 1,
                    "minor_version": 2
                },
                {
                    "type": "writers.gdal",
                    "filename": outtif,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }
    else: 
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
                },
                {
                    "type": "filters.dem",
                    "raster": dem_fp,
                    "limits": f"Z[{dem_low}:{dem_high}]"
                },
                {
                    "type": "filters.elm"
                },
                {
                    "type": "filters.outlier",\
                    "method": "statistical",\
                    "mean_k": mean_k,\
                    "multiplier": multiplier
                },
                # {
                #     "type": "filters.smrf",\
                #     "ignore": "Classification[7:7], NumberOfReturns[0:0], ReturnNumber[0:0]"
                # },
                {
                    "type": "filters.range",
                    "limits": "Classification[2:2]"
                },
                {
                    "type": "writers.las",
                    "filename": outlas,
                    "major_version": 1,
                    "minor_version": 2
                },
                {
                    "type": "writers.gdal",
                    "filename": outtif,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }

    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'dtm_pipeline'
    json_to_use = join(json_dir, f'{json_name}.json')

    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)

    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return outlas, outtif

def surface_models(laz_fp, outlas = '', outtif = '', user_dem = '', dem_low = 20, dem_high = 50, mean_k = 20, multiplier = 3, lidar_pc = 'yes'):
    """Use filters.dem, filters.mongo, filters.elm, filters.outlier, filters.smrf, and filters.range to filter the point cloud for surface models.

    Args:
        laz_fp (_type_): _description_
        outlas (str, optional): _description_. Defaults to ''.
        outtif (str, optional): _description_. Defaults to ''.
        user_dem (str, optional): _description_. Defaults to ''.
        dem_low (int, optional): _description_. Defaults to 20.
        dem_high (int, optional): _description_. Defaults to 50.
        mean_k (int, optional): _description_. Defaults to 20.
        multiplier (int, optional): _description_. Defaults to 3.
    
    Returns:
        _type_: Filepath to the terrain model.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    # os.chdir(in_dir)

    #create a filepath for the output las and tif file
    if outlas == '':
        outlas = join(in_dir, 'dsm.laz')
    if outtif == '':
        outtif = join(in_dir, 'dsm.tif')

    #set dem_fp
    dem_fp = join(in_dir, 'dem.tif')
    
    #download dem using download_dem() if user_dem is not provided
    if user_dem == '':
        dem_fp, crs, project = download_dem(laz_fp, dem_fp= dem_fp)
    else:
        shutil.copy(user_dem, dem_fp) #if user_dem is provided, copy the user_dem to dem_fp 

    if lidar_pc.lower() == 'yes':
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
                },
                {
                    "type": "filters.dem",
                    "raster": dem_fp,
                    "limits": f"Z[{dem_low}:{dem_high}]"
                },
                {
                    "type": "filters.mongo",\
                    "expression": {"$and": [\
                    {"ReturnNumber": {"$gt": 0}},\
                    {"NumberOfReturns": {"$gt": 0}} ] }
                },
                {
                    "type": "filters.elm"
                },
                {
                    "type": "filters.outlier",\
                    "method": "statistical",\
                    "mean_k": mean_k,\
                    "multiplier": multiplier
                },
                {"type": "filters.range",\
                "limits":"returnnumber[1:1]"
                },
                {
                    "type": "writers.las",
                    "filename": outlas,
                    "major_version": 1,
                    "minor_version": 2
                },
                {
                    "type": "writers.gdal",
                    "filename": outtif,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }
    else:
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
                },
                {
                    "type": "filters.dem",
                    "raster": dem_fp,
                    "limits": f"Z[{dem_low}:{dem_high}]"
                },
                {
                    "type": "filters.elm"
                },
                {
                    "type": "filters.outlier",\
                    "method": "statistical",\
                    "mean_k": mean_k,\
                    "multiplier": multiplier
                },
                {"type": "filters.range",\
                "limits":"Classification[2:6]"
                },
                {
                    "type": "writers.las",
                    "filename": outlas,
                    "major_version": 1,
                    "minor_version": 2
                },
                {
                    "type": "writers.gdal",
                    "filename": outtif,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }

    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'dsm_pipeline'
    json_to_use = join(json_dir, f'{json_name}.json')

    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
        
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return outlas, outtif