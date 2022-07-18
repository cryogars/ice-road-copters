# ice-road-copters
The goal of this program is to utilize existing USGS 3DEP high resolution topography data and Ames Stereo Pipeline (ASP) software to align snow-on airborne lidar points clouds to real world coordinates without the use of ground control points.
![heli_bsu](./docs/heli.png)

For our case, we used prior knowledge that HWY-21 running through our study site is kept snow-free for a majority of the year, thus making excellent virtual ground control points for post-processing in ASP.
![roads](./docs/roads.png)

###  Set up conda environment (Zach, I had to create an issue with ASP... conda was not creating an env for me... https://github.com/NeoGeographyToolkit/StereoPipeline/issues/372)
```
$ conda env create -f asp_3.1.0_linux_env.yaml
$ conda activate asp
$ conda install -c conda-forge pandas matplotlib geopandas python-pdal rasterio xarray whitebox
```
NOTE: This environment containing ASP will only work on linux operating systems.

### USGS 1m DEM tiles were downloaded along the HWY-21 using `wget` via https://apps.nationalmap.gov/downloader/