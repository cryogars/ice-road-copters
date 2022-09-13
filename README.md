# ice-road-copters :helicopter:

###  Setting things up! :hammer:
#### Downloading ASP precompiled binaries
1. Download latest stable build (Linux or OSx (3.1.0)): https://github.com/NeoGeographyToolkit/StereoPipeline/releases and unzip the folder into the ice-road-copters directory. There are different builds for each OS, but you may have to dig a little to find the OSx build, it gets updated less.
2. Rename this folder as `ASP` and remove the zipped file


#### Setting up Conda environment 

```
$ conda env create -f iceroad_env.yaml
$ conda activate iceroad
```

#### Running the code
from the ice-road-copters directory for example one can run:
```
$ python scripts/ice-road-pipeline.py <path-to-directory-of-laz-files> -e <path-to-user-supplied-reference-dem> -a <path-to-ASP-directory> -s <path-to-road-shapefile-to-clip-to>
```
NOTE: that this code assumes you are using a reference DEM (and heli lidar data) in ellipsoid. If your reference DEM is in Geoid you must set the `-g = True` flag.

#### Flags

```
-e user_dem      Path to user specifed DEM [one will be downloaded from py3dep if you don't supply one]
-d debug         turns on debugging logging [Options: True or False]
-a asp_dir       Directory with ASP binary files [Can be either ASP or ASP/bin directory]
-s shp_fp        Shapefile to align with [road shapefile to use to tie reference DEM to your point cloud]
-g geoid         Is the reference DEM in geoid? Will be auto set to True if you don't supply a DEM [Default: False]
```


###  Additional information :books:
The goal of this program is to utilize existing USGS 3DEP 1m topography data (via [py3dep](https://github.com/hyriver/py3dep)) and Ames Stereo Pipeline ([ASP](https://github.com/NeoGeographyToolkit/StereoPipeline)) software to accurately align snow-on airborne lidar point clouds to real world coordinates without the use of ground control points.

![heli_bsu](./docs/heli.png). We also provide an option for a user specified DEM (snow-off).

For our case, we used prior knowledge that HWY-21 running through our study site is kept snow-free for a majority of the year, thus making excellent virtual ground control points for post-processing in ASP.

![roads](./docs/roads.png)

After running pdal processing and ASP post-processing, we are able to generate accurate snow depth maps with sub cm road differences.

![snow](./docs/snow.jpeg)
