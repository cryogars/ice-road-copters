import os
from os.path import dirname, join
import json
import subprocess
import shutil
from snow_pc.common import download_dem, make_dirs

def return_filtering(laz_fp, out_fp = ''):
    """Use filters.mongo to filter out points with invalid returns.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.

    Returns:
        _type_: Filepath to the filtered point cloud file.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #create a filepath for the output las file
    if out_fp == '':
        out_fp = "returns_filtered.laz"

    #create a json pipeline for pdal
    json_pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": laz_fp
            },
            {
                "type": "filters.mongo",\
                "expression": {"$and": [\
                {"ReturnNumber": {"$gt": 0}},\
                {"NumberOfReturns": {"$gt": 0}} ] }
            },
            {
                "type": "writers.las",
                "filename": out_fp
            }
        ]
    }
    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'return_filtering'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp

def dem_filtering(laz_fp, user_dem = '', dem_low = 20, dem_high = 50, out_fp = ''):
    """Use filters.dem to filter the point cloud to the DEM. 

    Args:
        laz_fp (_type_): Filepath to the point cloud file.
        user_dem (str, optional): Filepath to the DEM file. Defaults to ''.
        dem_low (int, optional): Lower limit of the DEM. Defaults to 20.
        dem_high (int, optional): Upper limit of the DEM. Defaults to 50.

    Returns:
        _type_: Filepath to the filtered point cloud file.
    """

    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #set dem_fp
    dem_fp = join(in_dir, 'dem.tif')
    
    #download dem using download_dem() if user_dem is not provided
    if user_dem == '':
        dem_fp, crs, project = download_dem(laz_fp, dem_fp= dem_fp)
    else:
        shutil.copy(user_dem, dem_fp) #if user_dem is provided, copy the user_dem to dem_fp

    #create a filepath for the output las file
    if out_fp == '':
        out_fp = "dem_filtered.laz"


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
                "type": "writers.las",
                "filename": out_fp
            }
        ]
    }
    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'dem_filtering'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp

def elm_filtering(laz_fp, out_fp = ''):
    """Use filters.elm to filter the point cloud.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.

    Returns:
        _type_: Filepath to the filtered point cloud file.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #create a filepath for the output las file
    if out_fp == '':
        out_fp = "elm_filtered.laz"
    
    #create a json pipeline for pdal
    json_pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": laz_fp
            },
            {
                "type": "filters.elm"
            },
            {
                "type": "writers.las",
                "filename": out_fp
            }
        ]
    }
    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'elm_filtering'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp

def outlier_filtering(laz_fp, mean_k = 20, multiplier = 3, out_fp = ''):
    """Use filters.outlier to filter the point cloud.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.
        mean_k (int, optional): _description_. Defaults to 20.
        multiplier (int, optional): _description_. Defaults to 3.

    Returns:
        _type_: Filepath to the filtered point cloud file.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #create a filepath for the output las file
    if out_fp == '':
        out_fp = "outlier_filtered.laz"

    #create a json pipeline for pdal
    json_pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": laz_fp
            },
            {
                "type": "filters.outlier",\
                "method": "statistical",\
                "mean_k": mean_k,\
                "multiplier": multiplier
            },
            {
                "type": "writers.las",
                "filename": out_fp
            }
        ]
    }
    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'outlier_filtering'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp

def ground_segmentation(laz_fp, out_fp = '', out_fp2 = '', lidar_pc = 'yes'):
    """Use filters.smrf and filters.range to segment ground points.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.

    Returns:
        _type_: Filepath to the segmented point cloud file.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #create a filepath for the output laz and tif file
    if out_fp == '':
        out_fp = "ground_segmented.laz"
    if out_fp2 == '':
        out_fp2 = "ground_segmented.tif"

    if lidar_pc.lower() == 'yes':
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
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
                    "filename": out_fp
                },
                {
                    "type": "writers.gdal",
                    "filename": out_fp2,
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
                    "type": "filters.range",
                    "limits": "Classification[2:2]"
                },
                {
                    "type": "writers.las",
                    "filename": out_fp,
                    "major_version": 1,
                    "minor_version": 4
                },
                {
                    "type": "writers.gdal",
                    "filename": out_fp2,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }
    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'ground_segmentation'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp, out_fp2

def surface_segmentation(laz_fp, out_fp = '', out_fp2 = '', lidar_pc = 'yes'):
    """Use filters.range to segment the surface points.

    Args:
        laz_fp (_type_): Filepath to the point cloud file.

    Returns:
        _type_: Filepath to the segmented point cloud file.
    """
    #set the working directory
    in_dir = os.path.dirname(laz_fp)
    os.chdir(in_dir)

    #create a filepath for the output laz and tif file
    if out_fp == '':
        out_fp = "surface_segmented.laz"
    if out_fp2 == '':
        out_fp2 = "surface_segmented.tif"

    if lidar_pc.lower() == 'yes':    
        #create a json pipeline for pdal
        json_pipeline = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": laz_fp
                },
                {
                    "type": "filters.range",
                    "limits": "returnnumber[1:1]"
                },
                {
                    "type": "writers.las",
                    "filename": out_fp
                },
                {
                    "type": "writers.gdal",
                    "filename": out_fp2,
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
                {"type": "filters.range",\
                "limits":"Classification[2:6]"
                },
                {
                    "type": "writers.las",
                    "filename": out_fp,
                    "major_version": 1,
                    "minor_version": 4
                },
                {
                    "type": "writers.gdal",
                    "filename": out_fp2,
                    "resolution": 1.0,
                    "output_type": "idw"
                }
            ]
        }

    #create a directory to save the json pipeline
    json_dir =  join(in_dir, 'jsons')
    os.makedirs(json_dir, exist_ok= True)
    json_name = 'surface_segmentation'
    json_to_use = join(json_dir, f'{json_name}.json')
    #write json pipeline to file
    with open(json_to_use, 'w') as f:
        json.dump(json_pipeline, f)
    #run the json pipeline
    subprocess.run(["pdal", "pipeline", json_to_use])

    return out_fp, out_fp2