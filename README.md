# ice-road-copters
The goal of this program is to utilize existing USGS 3DEP high resolution topography data and Ames Stereo Pipeline (ASP) software to accurately align snow-on airborne lidar point clouds to real world coordinates without the use of ground control points.
![heli_bsu](./docs/heli.png)

For our case, we used prior knowledge that HWY-21 running through our study site is kept snow-free for a majority of the year, thus making excellent virtual ground control points for post-processing in ASP.
![roads](./docs/roads.png)


###  Setting up things up

#### Downloading ASP precompiled binaries
1. Download latest stable build (Linux or OSx): https://github.com/NeoGeographyToolkit/StereoPipeline/releases and unzip the folder into the ice-road-copters directory
2. Rename this folder as `ASP` and remove the zipped file


#### Setting up Conda environment 

```
$ conda env create -f iceroad_env.yaml
$ conda activate iceroad
```


### USGS 1m DEM tiles were downloaded along the HWY-21 using `wget` via https://apps.nationalmap.gov/downloader/