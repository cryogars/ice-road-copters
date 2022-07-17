# ice-road-copters
The goal of this program is to utilize existing USGS 3DEP high resolution topography data and Ames Stereo Pipeline (ASP) software to align snow-on airborne lidar points clouds to real world coordinates without the use of ground control points.
![heli_bsu](./docs/heli.png)

###  Set up conda environments
#### Lidar processing env
```
$ conda env create -f environment.yml
$ conda activate heli
```


#### ASP env
```
$ conda env create -f asp_env.yml
$ conda activate asp
$ which stereo
```



