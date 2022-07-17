# quick steps for us to get started...

### 0. Get the raw lidar data
Heli data: `SNOWDATA/Nah/2022_IDALS`


###  1. Set up a conda environment
```
$ conda env create -f environment.yml
$ conda activate heli
```


###  2. Set up asp conda environment
```
$ conda env create -f asp_env.yml
$ conda activate asp
$ which stereo
```



### 3. let's do it! 

Dom's code: https://github.com/granitehills14/Airborne-Laser-Scanning


ASP example wrappers: https://github.com/friedrichknuth/hsfm/blob/master/hsfm/asp/asp.py
